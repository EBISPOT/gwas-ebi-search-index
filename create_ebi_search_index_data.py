import cx_Oracle
import contextlib
import sys
import os
import argparse
import subprocess
import json
from subprocess import Popen, PIPE
from datetime import date
from gwas_db_connect import DBConnection


class ebiSearchIndexData(object):

    getTotalStudiesSQL = '''
        SELECT COUNT(DISTINCT (S.ACCESSION_ID)) 
        FROM STUDY S, HOUSEKEEPING H, PUBLICATION P 
        WHERE S.HOUSEKEEPING_ID=H.ID and H.IS_PUBLISHED=1 
            and S.PUBLICATION_ID=P.ID
    '''

    getDataSQL = '''
        SELECT P.PUBMED_ID, S.ID, S.ACCESSION_ID, S.INITIAL_SAMPLE_SIZE, AU.LAST_NAME, AU.INITIALS, 
            DT.TRAIT AS REPORTED_TRAIT, listagg(ET.SHORT_FORM, ', ') WITHIN GROUP (ORDER BY ET.SHORT_FORM) 
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

    findMissingStudiesSQL = '''
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

    def __init__(self, connection, database):
        # self.outputDir = outputDir
        self.database = database

        try:
            with contextlib.closing(connection.cursor()) as cursor:
                ######################
                # Get all study data
                ######################
                cursor.execute(self.getDataSQL)
                data = cursor.fetchall()
                # print("Num of Studies: ", len(study_data))
                
                # Get first X rows
                self.data = data[:2]
                
                # Get row at index 99 in list
                # self.data = data[99:100]
                
                # Get all data
                # self.data = data


                #####################
                # Get total studies
                #####################
                cursor.execute(self.getTotalStudiesSQL)
                total_studies = cursor.fetchone()[0]
                # print("Total studies: ", total_studies)
                self.total_studies = total_studies

        except(cx_Oracle.DatabaseError, exception):
            print(exception)


    def data_check(self):
        """ Confirm that total number of results from getDataSQL query match getTotalStudies query. 
        If not equal, there are studies that are missing reported or mapped trait annotations.
        """
        if self.total_studies != len(self.data):
            print('[Warning] Number of studies do not match.')
            self.__find_data_errors()


    def __find_data_errors(self):
        """ Find studies that are missing reported trait or mapped trait annotations."""
        try:
            with contextlib.closing(connection.cursor()) as cursor:
                #####################################
                # Find studies missing annotations
                ####################################
                cursor.execute(self.findMissingStudiesSQL)
                incorrect_studies = cursor.fetchall()
                print(incorrect_studies)
                # QUESTION: Should a mismatch in study annotations cause Data Release or this script to Fail? The reported trait is 
                # used in the name of an entry in the EBI Search Index snippet
                # TODO: Send email to gwas-curators with information about incorrect_studies
        except(cx_Oracle.DatabaseError, exception):
            print(exception)


    def formatData(self):
        """ Format data into EBI Search Index JSON format. 

        Returns:
            JSON file
        """
        
        header = self.__create_data_header()
        print("** Header: ", header)

        # dataObj = {"entries": [{"fields": [], "cross_references": []}]}
        entries_list = []
        # entries_dict = {"entries": ""}

        # print(self.data)
        for pubmed_Id, study_Id, accession_Id, initial_sample_size, last_name, initials, reported_trait, efo_short_form in self.data:

        # for row in self.data:
            # print("Row: ", row[7:])
            # print("\nRow: ", row)

            
            entry = {"fields": "", "cross_references": ""}
            fields_list = []
            id_field = {"name": "id", "value": ""}
            name_field = {"name": "name", "value": ""}
            description_field = {"name": "description", "value": ""}
            cross_references_list = []
            cross_references_dict = {}

            # Populate individual "fields" dictionaries
            id_field['value'] = accession_Id
            print("** ID Field: ", id_field)

            name_field['value'] = 'GWAS: ' + reported_trait + ' (' + accession_Id + ')'


            description_field['value'] = 'Study of ' + initial_sample_size + ' published by ' + last_name + ' ' + initials + ', PMID: ' + pubmed_Id


            # Add "id_field" dictionary to "fields_list"
            fields_list.append(id_field)
            fields_list.append(name_field)
            fields_list.append(description_field)
            print("** FL: ", fields_list)

            # Add "fields_list" list to "entry" dictionary
            entry['fields'] = fields_list
            print("** Entry: ", entry)


            # Add entry dictionary to "entries_list"
            entries_list.append(entry)


        print("\n\n** All Entry List: ", entries_list)

        # entries_dict['entries'] = entries_list
        # print("\n\n** All Entries: ", entries_dict)

        header['entries'] = entries_list
        print("\n\n** All Data: ", header)

        print("\n\nAll Data JSON: ", json.dumps(header))


    def __create_data_header(self):
        """Create data file header.

        Returns:
            dict: Dictionary containing header data.

        Example: 
            {"entry_count": 9143, 
            "name": "GWAS Catalog",
            "release": "22-04-2020",
            "release_date": "22-04-2020"}
        """
        header = {"entry_count": len(self.data), "name": "GWAS Catalog", "release": "", "release_date": "", "entries": ""}
                
        date_today = date.today().strftime('%d-%m-%Y')
        header['release'] = date_today
        header['release_date'] = date_today

        return header



if __name__ == '__main__':

    # Parsing command line arguments:
    parser = argparse.ArgumentParser()
    parser.add_argument('--releaseDB', type=str, help='Name of the database for extracting study data.')
    # parser.add_argument('--outputDir', type=str, help='Path to data directory.')
    # parser.add_argument('--emailRecipient', type=str, help='Email address where the notification is sent.')
    args = parser.parse_args()

    database = args.releaseDB
    # outputDir = args.outputDir
    # emailRecipient = args.emailRecipient

    # Check if output directory exists:
    # if not os.path.isdir(outputDir):
    #     print('[Error] No valid staging directory provided. Exiting.')
    #     sys.exit(1)


    # Open connection:
    db_object = DBConnection.gwasCatalogDbConnector(database)
    connection = db_object.connection

    # Get published studies from database
    ebiSearchIndexDataObj = ebiSearchIndexData(connection, database)

    # Check data integrity
    ebiSearchIndexDataObj.data_check()

    # Format data
    ebiSearchIndexDataObj.formatData()



