#!/usr/bin/env python

"""Create results report."""

import argparse
from aplanat import report
import pandas as pd
import os
import numpy as np
from bokeh.resources import INLINE

def read_patient_info(file, barcode):
    """
    Read patient info file and return relevant row for the barcode
    
    Args:
        file (str): path to the patient info file
        barcode (str): barcode identifier for the patient sample
        
    Returns:
        list: list of pandas Series with patient info
    """
    #funtion parsing CSV patient file and looking up info for relevant barcode
    info=pd.read_excel(file, usecols=range(0,11))
    #return row with a barcode as a Series
    if barcode=="discontinued":
        relevant_row=[]
        relevant_rows=info.loc[info['Status'] == barcode]
        for index,row in relevant_rows.iterrows():
            single_row=row.transpose()
            single_df=single_row.reset_index()
            single_df.columns=['Metadata', 'Sample Information']
            relevant_row.append(single_df)
    else:
        relevant_rows=info.loc[info['Barcode'] == barcode]
        #move row names into a column
        relevant_row=relevant_rows.transpose()
        relevant_row.index.name = 'Metadata'
        relevant_row.reset_index(inplace=True)
        #rename columns
        relevant_row.columns=['Metadata', 'Sample Information']
    return relevant_row

def read_abundance_results(file):
    """
    Read abundance results csv file and return top 3 results

    Args:
        file (path): path to the abundance results file    

    Returns:
        pandas.DataFrame: top 3 abundance results
    """
    # The abundance results file is a CSV file
    abundance_results = pd.read_csv(file)
    
    # Rename the columns
    abundance_results.columns = ['Detected Species', 'Relative Abundance (%)', 'Number of Reads']
    
    # Round the abundance results
    abundance_results['Relative Abundance (%)'] = abundance_results['Relative Abundance (%)'].apply(np.around)
    
    # Return the top 3 abundance results
    abundance_results_t3 = abundance_results.head(n=3)
    
    return abundance_results_t3

def process_controls(ctrl):
    """
    Process controls files

    Args:
        control (path): Path to the positive or negative control file

    Returns:
        dataframe: Pandas dataframe with the control results
    """
    ctrl = ctrl[1:-1]
    if ctrl == '' or ctrl == 'None':
        control = None
    else:
        control = read_abundance_results(ctrl)
            

    return control

def parse_args():
    """Run the entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--infile", default='unknown',
        help='Table file with classification and abundance results')
    parser.add_argument(
        "--output", default='unknown',
        help="Report output file name")
    parser.add_argument(
        "--barcode", default='discontinued',
        help="barcode identifier for the patient sample")
    parser.add_argument(
        "--info", required=True,
        help="Experiment information file mapping patient ID with a barcode")
    parser.add_argument(
        "--revision", default='unknown',
        help="git branch/tag of the executed workflow")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit of the executed workflow")
    parser.add_argument(
        "--demux", default='unknown',
        help='demultiplexing method and software version')
    parser.add_argument(
        "--clustering_size", default='unknown',
        help="Amount of reads used for UMAP HDBSCAN clustering")
    parser.add_argument(
        "--positive", default='unknown',
        help="File path for positive control results")
    parser.add_argument(
        "--negative", default='unknown',
        help="File path for negative control results")
    parser.add_argument(
        "--reads_count", default='0',
        help="Reads count after quality control")
    parser.add_argument(
        "--kit", default='unknown',
        help="Kit used for barcoding and demultiplexing")
    parser.add_argument(
        "--report_template",
        help="path to report template")
    parser.add_argument(
        "--logo",
        help="custom logo")
    parser.add_argument(
        "--run_id", default='unknown',
        help="Run ID of the sequencing run")
    parser.add_argument(
        "--seq_start", default='unknown',
        help="Start time of the sequencing run")
    
    args = parser.parse_args()

    return(args)

def main(args):
    
    if args.barcode=="discontinued":

        metadata_table_list=read_patient_info(args.info, args.barcode)
                
        for patient in metadata_table_list:
            # Restructure the metadata table
            restructured=[]
            for index, row in patient.iloc[[0,2,1,4,6,7,10]].iterrows():
                restructured.append(": ".join([str(row['Metadata']), str(row['Sample Information'])]))
            restructured.insert(4, " ".join(["Sequencing start:", args.seq_start]))
            rest_df=pd.DataFrame(list(zip(restructured[:4],restructured[4:])), columns=['Sample Information', 'Time Stamps'])

            # Generate the report
            title="Patient " + patient.iloc[0,1] + " Report"
            reprt = report.UoSReport(
                title=title, report_template=args.report_template, about=False, style='UoS', logo=args.logo)

            section=reprt.add_section()
            section.markdown('''
            ### Sample Information
            ''')

            section.table(rest_df)

            section=reprt.add_section()

            assay_type=patient.iloc[4]['Sample Information']
            if assay_type == '16S':
                assay_info = 'Bacterial 16s'
            else:
                assay_info = 'Fungal ITS2'

            section.markdown('''
            <br/>
            ### Results

            <font color="red">**{0} rRNA NOT detected**</font>
            '''.format(assay_info))

            reprt.write("patient_report_" + str(patient.iloc[1,1]) + ".html")
    else:        
        metadata_table=read_patient_info(args.info, args.barcode)
        restructured=[]
        for index, row in metadata_table.iloc[[0,2,1,4,6,7,10]].iterrows():
            restructured.append(": ".join([str(row['Metadata']), str(row['Sample Information'])]))
        restructured.insert(4, " ".join(["Sequencing start:", args.seq_start]))
        rest_df=pd.DataFrame(list(zip(restructured[:4],restructured[4:])), columns=['Sample Information', 'Time Stamps'])

        # Create the title for the report
        title="Patient " + metadata_table.iloc[0,1] + " Report"

        if args.infile == "input.1":
            results_table = None
        else:
            results_table=read_abundance_results(args.infile)

        positive=process_controls(args.positive)
        negative=process_controls(args.negative)

        reprt = report.UoSReport(
            title=title, workflow="NanoCLUST", report_template=args.report_template,
            revision=args.revision, commit=args.commit, style='UoS', logo=args.logo)

        section=reprt.add_section()
        section.markdown('''
        ### Sample Information
        ''')

        section.table(rest_df)

        section=reprt.add_section()

        section.markdown('''
        <br/>
        ### Results

        Total reads in this sample: {0}
        '''.format(args.reads_count))

        assay_type=metadata_table.loc[metadata_table['Metadata'] == 'Assay', 'Sample Information'].iloc[0]

        if assay_type == '16S':
            database_info="16s bacterial sequencing results were compared against 16S & 18S database, build 18 Jan 2022."
            infection_type='bacterial'
            assay_info='Bacterial 16s'
        else:
            database_info="ITS2 fungal sequencing results were compared against ITS2 database, build 15 Mar 2022."
            infection_type='fungal'
            assay_info='Fungal ITS2'

        if results_table is not None:
            section.table(results_table, classes=['highlighted', 'larger-first-column'])
        else:
            section.markdown('''
            <font color="red">**{0} rRNA NOT detected**</font>
            '''.format(assay_info))
        

        section=reprt.add_section()
        section.markdown('''
        <br/>
        ### Run QC


        **NEGATIVE CONTROL**
        ''')
        comment=0
        if negative is not None:
            comment+=10
            section.markdown('''
            Total reads in negative control: {0} 
            '''.format(negative['Number of Reads'].sum()))

            section.table(negative, classes='larger-first-column')
        else:
            section.markdown('''
            No species detected in negative control.
            ''')

        section.markdown('''
        <br/>
        **POSITIVE CONTROL**
        ''')

        if positive is not None:
            comment+=1
            section.markdown('''
            Total reads in positive control: {0}
            '''.format(positive['Number of Reads'].sum()))
            
            section.table(positive, classes='larger-first-column')
        else:
            comment+=2
            section.markdown('''
            No species detected in positive control.
            ''')

        if comment == 1:
            section.markdown('''
            comment: <font color="green">QC for this sample was **successful**</font>
            ''')
        elif comment == 11:
            section.markdown('''
            comment: <font color="orange">QC for this sample shows the presence of {0} reads in the negative control. Check if there is overlap with any detected pathogen in the sample that may indicate contamination.</font>
            '''.format(infection_type))
        else:
            section.markdown('''
            comment: <font color="red">QC for this sample has **failed** and results cannot be validated</font>
            ''')

        section=reprt.add_section()
        run_id=args.run_id
        barcoding_kit=args.kit
        print(barcoding_kit)
        demux_method=args.demux
        species_database=metadata_table['Sample Information'].iloc[4]
        clustering_size=args.clustering_size

        run_params=pd.DataFrame(list(zip(['Run ID: '+str(run_id), 'Barcoding kit: '+str(barcoding_kit), 'Demultiplex method: '+str(demux_method)], ['Species database: '+str(species_database), 'Clustering size: '+str(clustering_size), 'Sample barcode: '+str(metadata_table.iloc[1,1])])), columns=['GridIon properties', 'NanoCLUST properties'])
        #run_params.reset_index(drop=True, inplace=True)

        section.markdown('''
        ### Run parameters
        ''')
        section.table(run_params)
        section.markdown('''
        Sample was sequenced on a ONT GridION Mk1. 
        Sequencing data was processed and analysed using a custom nanoclust pipeline.
        {6}

        '''.format(run_id, barcoding_kit, demux_method, species_database, clustering_size, metadata_table.iloc[1,1], database_info))

        #write report
        reprt.write(args.output + "_" + str(metadata_table.iloc[1,1]) + ".html")


if __name__ == "__main__":
    args = parse_args()

    main(args)