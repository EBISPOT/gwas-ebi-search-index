import cx_Oracle
import contextlib
import sys
import os
import argparse
import subprocess
import json
from tqdm import tqdm
from subprocess import Popen, PIPE
import smtplib
from os.path import basename
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.application import MIMEApplication
from datetime import datetime, date
from gwas_db_connect import DBConnection


class EbiSearchIndexData:

    TOTAL_STUDIES_SQL = '''
        SELECT COUNT(DISTINCT (S.ACCESSION_ID)) 
        FROM STUDY S, HOUSEKEEPING H, PUBLICATION P 
        WHERE S.HOUSEKEEPING_ID=H.ID and H.IS_PUBLISHED=1 
            and S.PUBLICATION_ID=P.ID
    '''

    ALL_STUDY_DATA_SQL = '''
        SELECT P.PUBMED_ID, P.TITLE, P.PUBLICATION_DATE, P.PUBLICATION, S.ID, S.ACCESSION_ID, S.INITIAL_SAMPLE_SIZE, AU.FULLNAME, 
            DT.TRAIT AS REPORTED_TRAIT, listagg(ET.SHORT_FORM, ', ') WITHIN GROUP (ORDER BY ET.SHORT_FORM) 
        FROM STUDY S, HOUSEKEEPING H, AUTHOR AU, PUBLICATION P, STUDY_DISEASE_TRAIT SDT, DISEASE_TRAIT DT, STUDY_EFO_TRAIT SETR, EFO_TRAIT ET
        WHERE S.PUBLICATION_ID=P.ID and S.HOUSEKEEPING_ID=H.ID and H.IS_PUBLISHED=1 and P.FIRST_AUTHOR_ID=AU.ID 
            and S.ID=SDT.STUDY_ID and SDT.DISEASE_TRAIT_ID=DT.ID
            and  S.ID=SETR.STUDY_ID and SETR.EFO_TRAIT_ID=ET.ID 
            and S.PUBLICATION_ID IN (
                SELECT P.ID
                FROM PUBLICATION P
            )
        GROUP BY P.PUBMED_ID, P.TITLE, P.PUBLICATION_DATE, P.PUBLICATION, S.ID, S.ACCESSION_ID, S.INITIAL_SAMPLE_SIZE, AU.FULLNAME, DT.TRAIT
    '''

    MISSING_STUDIES_SQL = '''
        SELECT P.PUBMED_ID, S.ID AS STUDY_ID, S.ACCESSION_ID 
        FROM STUDY S, HOUSEKEEPING H, PUBLICATION P 
        WHERE S.HOUSEKEEPING_ID=H.ID and H.IS_PUBLISHED=1 
            and S.PUBLICATION_ID=P.ID
        MINUS 
        SELECT P.PUBMED_ID, S.ID AS STUDY_ID, S.ACCESSION_ID 
        FROM STUDY S, HOUSEKEEPING H, AUTHOR AU, PUBLICATION P, STUDY_DISEASE_TRAIT SDT, DISEASE_TRAIT DT, STUDY_EFO_TRAIT SETR, EFO_TRAIT ET
        WHERE S.PUBLICATION_ID=P.ID and S.HOUSEKEEPING_ID=H.ID and H.IS_PUBLISHED=1 and P.FIRST_AUTHOR_ID=AU.ID 
            and S.ID=SDT.STUDY_ID and SDT.DISEASE_TRAIT_ID=DT.ID 
            and  S.ID=SETR.STUDY_ID and SETR.EFO_TRAIT_ID=ET.ID 
            and S.PUBLICATION_ID IN (
                SELECT P.ID
                FROM PUBLICATION P
            )
        GROUP BY P.PUBMED_ID, S.ID, S.ACCESSION_ID, S.INITIAL_SAMPLE_SIZE, AU.LAST_NAME, AU.INITIALS, DT.TRAIT
    '''

    studies_missing_data = []

    def __init__(self, connection, database, output_dir, logs_dir, email_recipient):
        self.database = database
        self.output_dir = output_dir
        self.logs_dir = logs_dir
        self.email_recipient = email_recipient

        try:
            with contextlib.closing(connection.cursor()) as cursor:
                ######################
                # Get all study data
                ######################
                cursor.execute(self.ALL_STUDY_DATA_SQL)
                data = cursor.fetchall()
                self.data = data
                # self.data = data[:20]

                #####################
                # Get total studies
                #####################
                cursor.execute(self.TOTAL_STUDIES_SQL)
                total_studies = cursor.fetchone()[0]
                self.total_studies = total_studies

        except cx_Oracle.DatabaseError as exception:
            print(exception)


    def data_check(self, connection):
        ''' Confirm that total number of results from ALL_STUDY_DATA_SQL query match TOTAL_STUDIES_SQL query. 
        If not equal, there are studies that are missing reported or mapped trait annotations.
        '''
        if self.total_studies != len(self.data):
            print('[Warning] Number of studies do not match.')
            # self.studies_missing_data.append('[Warning] Number of studies do not match.')
            self._find_data_errors(connection)


    def _find_data_errors(self, connection):
        ''' Find studies that are missing reported trait or mapped trait annotations.'''
        try:
            with contextlib.closing(connection.cursor()) as cursor:
                #####################################
                # Find studies missing annotations
                ####################################
                cursor.execute(self.MISSING_STUDIES_SQL)
                incorrect_studies = cursor.fetchall()
               
                with open(self.logs_dir + 'missingStudies.txt', 'w') as file:
                    print(incorrect_studies, file=file)
                
                # TODO: Send email to gwas-curators with information about incorrect_studies
                # report_string += '\n'.join(map('\t-- {}'.format, data))
                accessions_missing_trait_annotations = [study[2] for study in incorrect_studies]

                # self.studies_missing_data.append(', '.join(['Accession: ' + acc for acc in accessions_missing_trait_annotations]))
                formatted_accession_list = ['-- Accession: ' + acc for acc in accessions_missing_trait_annotations]

                self.studies_missing_data.append('Number of studies do not match. Missing trait information for accessions: ')
                self.studies_missing_data.append(', '.join(accessions_missing_trait_annotations))

                # for accession in formatted_accession_list:
                #     self.studies_missing_data.append(accession)

        except cx_Oracle.DatabaseError as exception:
            print(exception)


    def format_data(self):
        ''' Format data into EBI Search Index JSON format. 

        Returns:
            JSON file
        '''
        
        data_file_obj = self._populate_header_data()

        entries_list = []

        for pmid, title, date, journal_name, study_Id, accession_Id, initial_sample_size, author_fullname, reported_trait, efo_short_form in tqdm(self.data):
            #########################
            # Define local variables
            #########################        
            entry = {'fields': '', 'cross_references': ''}
        
            fields_list = []
            id_field = {'name': 'id', 'value': ''}
            reported_trait_field = {'name': 'reported_trait', 'value': ''}
            description_field = {'name': 'description', 'value': ''}

            url_published_field = {'name': 'url', 'value': ''}
            url_published_field_prefix = 'https://www.ebi.ac.uk/gwas/studies/'
            # url_prepublished_field_prefix = 'ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/'

            first_author_field = {'name': 'first_author', 'value': ''}
            publication_title_field = {'name': 'publication_title', 'value': ''}
            journal_name_field = {'name': 'journal_name', 'value': ''}
            publication_date_field = {'name': 'publication_date', 'value': ''}
            
            cross_references_list = []
            pmid_xref = {'dbkey': '', 'dbname': 'PUBMED'}
            efo_xref = {'dbkey': '', 'dbname': 'EFO'}


            #############################################
            # Populate individual 'fields' dictionaries
            #############################################
            if accession_Id is None:
                self.studies_missing_data.append('PMID: ' + pmid + ' has a study missing an accession identifier.')
                continue
            id_field['value'] = accession_Id
            url_published_field['value'] = url_published_field_prefix + accession_Id


            if reported_trait is None:
                self.studies_missing_data.append('Accession: ' + accession_Id + ' is missing a reported trait annotation.')
                continue
            reported_trait_field['value'] = reported_trait


            if initial_sample_size is None:
                self.studies_missing_data.append('Accession: ' + accession_Id + ' is missing inital sample size information.')
                description_field['value'] = 'Study published by ' + author_fullname + ', PMID: ' + pmid
            else:
                description_field['value'] = 'Study of ' + initial_sample_size


            if author_fullname is None:
                self.studies_missing_data.append('Accession: ' + accession_Id + ' is missing the First Author name.')
                continue
            first_author_field['value'] = author_fullname


            if title is None:
                self.studies_missing_data.append('Accession: ' + accession_Id + ' is missing a publication title.')
                continue
            publication_title_field['value'] = title


            if journal_name is None:
                self.studies_missing_data.append('Accession: ' + accession_Id + ' is missing a journal name.')
                continue
            journal_name_field['value'] = journal_name


            if date is None:
                self.studies_missing_data.append('Accession: ' + accession_Id + ' is missing a publication date.')
                continue
            publication_date_field['value'] = date.strftime('%Y/%m/%d')


            # Add to 'fields_list'
            fields_list.append(id_field)
            fields_list.append(reported_trait_field)
            fields_list.append(description_field)
            fields_list.append(first_author_field)
            fields_list.append(publication_title_field)
            fields_list.append(journal_name_field)
            fields_list.append(publication_date_field)
            fields_list.append(url_published_field)

            # Add 'fields_list' list to 'entry' dictionary
            entry['fields'] = fields_list


            ######################################################
            # Populate individual 'cross_references' dictionaries
            ######################################################
            if pmid is None:
                pmid_xref['dbkey'] = ''
            else:
                pmid_xref['dbkey'] = pmid
                # Add xrefs to 'cross_references' list
                cross_references_list.append(pmid_xref)

            if efo_xref is None:
                efo_xref['dbkey'] = ''
            else:
                # Account for studies that have more than one EFO trait assigned
                for i in range((len(efo_short_form.split(', ')))):
                    efo_xref['dbkey'] = efo_short_form.split(', ')[i]
                    
                    # Add xrefs to 'cross_references' list
                    cross_references_list.append(efo_xref)

                    # Clear dict
                    efo_xref = {'dbkey': '', 'dbname': 'EFO'}


            entry['cross_references'] = cross_references_list


            # Add entry dictionary to 'entries_list'
            entries_list.append(entry)


        # Add all entries to data file object 
        data_file_obj['entries'] = entries_list


        ###################
        # Save data file
        ###################
        with open(self.output_dir + 'studies.json', 'w') as file:
            json.dump(data_file_obj, file)


        ###################
        # Save log file
        ###################
        # TODO: Decide whether to send errors as file attachment or body of email 
        with open(self.logs_dir + 'logs.txt', 'w') as log_file:
            print(self.studies_missing_data, file=log_file)

        # Send email of errors as body of email
        return self.studies_missing_data


    def _populate_header_data(self):
        '''Add header data to data_file_obj.

        Returns:
            dict: Dictionary containing header data.

        Example: 
            {'entry_count': 9143, 
            'name': 'GWAS Catalog',
            'release': '22-04-2020',
            'release_date': '22-04-2020'}
        '''
        header = {'entry_count': len(self.data), 'name': 'GWAS Catalog', 'release': '', 'release_date': '', 'entries': ''}
                
        date_today = self._get_timestamp()
        header['release'] = date_today
        header['release_date'] = date_today

        return header


    def _get_timestamp(self):
        ''' 
        Get timestamp of current date. 
        '''
        return date.today().strftime('%d-%m-%Y')


    def send_email_report(self, studies_with_errors, email_addresses):
        ''' Send email with information about studies missing curation information '''
        try:
            mailBody = 'Subject: Data Release report - Published studies missing data annotations\nTo: {}\n{}'.format(email_addresses, studies_with_errors)
            p = Popen(["/usr/sbin/sendmail", "-t", "-oi", self.email_recipient], stdin=PIPE)
            p.communicate(mailBody.encode('utf-8'))
        except OSError as e:
            print(e) 


    def format_error_data(self, data):
        ''' Format data to add in body of email '''

        # Today's date
        date_today = self._get_timestamp()

        report_string = ('Report of publications or studies with missing curation information\n\n'
                            '[Info] Date of run: {}\n'
                            '[Info] Source database: {}\n\n'
                            '[Info] Errors found with:\n'
                            ''.format(date_today, self.database))

    
        if len(data) > 0:
            report_string += '\n'.join(map('\t-- {}'.format, data))
        else:
            report_string += ''.join('\nNo data issues found.')

        report_string += '\n\nThis report was created by the "create_ebi_search_index_data.py" script run by the "Data Release(prod)" Bamboo build plan.'

        return report_string

    

    def send_email_report_attachment(self, report_filename, email_recipient):
        ''' Email error report file '''
        
        # Today's date
        date_today = self._get_timestamp()

        with open(self.logs_dir + report_filename, "rb") as file:
            part = MIMEApplication(
                file.read(),
                Name=basename(report_filename)
            )
        
        # create a text/plain message
        msg = MIMEMultipart()

        # After the file is closed
        part['Content-Disposition'] = 'attachment; filename="%s"' % basename(report_filename)
        msg.attach(part)


        # create headers
        me = 'email_recipient'
        # you = ['gwas-dev-logs@ebi.ac.uk', 'gwas-curator@ebi.ac.uk']
        you = [email_recipient]
        msg['Subject'] = 'Published studies missing data ' +  date_today
        msg['From'] = me
        msg['To'] = ", ".join(you)

        # send the message via our own SMTP server, but don't include the envelope header
        s = smtplib.SMTP('localhost')
        s.sendmail(me, you, msg.as_string())
        s.quit()


def main():
    # Parsing command line arguments:
    parser = argparse.ArgumentParser()
    parser.add_argument('--release_db', type=str, help='Name of the database for extracting study data.')
    parser.add_argument('--output_dir', type=str, help='Path to data directory.')
    parser.add_argument('--logs_dir', type=str, help='Path to logs directory.')
    parser.add_argument('--email_recipient', type=str, help='Email address where the notification is sent.')
    args = parser.parse_args()

    database = args.release_db
    output_dir = args.output_dir
    logs_dir = args.logs_dir
    email_recipient = args.email_recipient

    # Check if output directory exists:
    if not os.path.isdir(output_dir):
        print('[Error] No valid data directory provided. Exiting.')
        sys.exit(1)

    # Check if logs directory exists:
    if not os.path.isdir(logs_dir):
        print('[Error] No valid logs directory provided. Exiting.')
        sys.exit(1)


    # Open connection:
    db_object = DBConnection.gwasCatalogDbConnector(database)
    connection = db_object.connection

    # Get published studies from database
    ebi_search_index_data_obj = EbiSearchIndexData(connection, database, output_dir, logs_dir, email_recipient)

    # Check data integrity
    ebi_search_index_data_obj.data_check(connection)

    # Format data
    error_data = ebi_search_index_data_obj.format_data()

    # Format error data
    formatted_report_data = ebi_search_index_data_obj.format_error_data(error_data)

    # Email error data
    ebi_search_index_data_obj.send_email_report(formatted_report_data, email_recipient)

    # Email report of missing curation information for studies as an attachment
    # ebi_search_index_data_obj.send_email_report_attachment('logs.txt', email_recipient)


if __name__ == '__main__':
    main()