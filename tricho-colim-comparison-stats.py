# convert to csv file
# extract proteins based on keywords in protein group file
# find peptides for these proteins in the psms file
# return table of intensities, sequence, protein, and sample for each peptide
# group and average the samples, with standard deviations
# plot the averages for each peptide

def convertCSV(directory):
    ''' converts .txt datafiles to .csv for easier reading into pandas
    '''
    import csv
    import os
    for dataFile in os.listdir(directory):
        if dataFile.endswith(".txt"):
            dataFile = os.path.join(directory, dataFile)
            csvFile = str(dataFile)[:-4] + ".csv"
        
            inTxt = csv.reader(open(dataFile, "rb"), delimiter = '\t')
            outCsv = csv.writer(open(csvFile, 'wb'))
            outCsv.writerows(inTxt)
        
def getProteinsFromKeyword(directory, keywords, description):
    ''' extracts protein descriptions that match keywords in the keyword list, and writes them to file
    '''
    import os
    import pandas as pd
    import csv
        
    outFile = os.path.join(directory, description + 'proteinList.csv')
    outCsv = csv.writer(open(outFile, 'wb'))
    protList = []
    
    for dataFile in os.listdir(directory):
        if dataFile.endswith("proteingroups.csv"):
            dataFile = os.path.join(directory, dataFile)
            data = pd.read_csv(dataFile)
            data = data.set_index('Accession')
            for keyword in keywords:
                for entry in data.index:
                    protein = data.at[entry,'Description']
                    if keyword in protein and protein not in protList:
                        outCsv.writerow([entry, protein])
                        protList.append(protein)
                
def getPeptides(directory, description):
    ''' extracts psm information from protein accessions extracted from getProteinsFromKeyword
        writes psm information to file
    '''
    import os
    import pandas as pd
    import csv
    
    outFile = os.path.join(directory, description + 'matchingpsms.csv')
    outCsv = csv.writer(open(outFile, 'wb'))
    
    protListFile = os.path.join(directory, description +  "proteinList.csv")
    accessions = pd.read_csv(protListFile, index_col=0, header=None)
    accessions = accessions.index.tolist()
    
    for dataFile in os.listdir(directory):
        if dataFile.endswith("_psms.csv"):
            dataFile = os.path.join(directory, dataFile) 
            data = pd.read_csv(dataFile)
            indeces = data.index.tolist()
            for entry in indeces:
                try:
                    accession = data.at[entry, 'Protein Group Accessions']
                    if int(accession) in accessions:
                        intensity = data.at[entry, 'Intensity']
                        sequence = data.at[entry, 'Sequence']
                        description = data.at[entry, 'Protein Descriptions']
                        sample = data.at[entry, 'Spectrum File']
                        sample = dataFile[-10]
                        modifications = data.at[entry, 'Modifications']
                        outCsv.writerow([accession, description, sequence, intensity, modifications, sample])
                except ValueError:
                    continue

def makeComparisons (directory, description):
    import pandas as pd
    import os
    import csv
    import numpy as np
    import scipy.stats as stats
    
    # read the psms file and group by sequence
    psmsListFile = os.path.join(directory, description + "matchingpsms.csv")
    psms = pd.read_csv(psmsListFile, names = ["Accession", "Description", "Sequence", "Intensity", "Modifications", "Sample"])
    groupedBySequence = psms.groupby(["Sequence"]) # first index in multiindex
    
    # setting up the write file
    outFile = os.path.join(directory, description + 'significantDifferences.csv')
    outCsv = csv.writer(open(outFile, 'wb'))  
    outCsv.writerow(["Accession", "Description", "Sequence", "Sample Group", "Mean", "Stdev", "colimited_mean/mean"])
    
    outFile2 = os.path.join(directory, description + "anova.csv")
    outCsv2 = csv.writer(open(outFile2, 'wb'))
    outCsv2.writerow(["Accession", "Sequence", "anovaf", "anovap", "toControl", "toP", "toFe"])
   
    for sequence, new_df in groupedBySequence:
        accession = new_df.iloc[0]["Accession"]
        description = new_df.iloc[0]["Description"]
        
        # calculate anova, select only sequences with statistically significant anova and colimitation t test values
        forStats = new_df.groupby("Sample")
        if len(forStats) >3:
            samples = []
            for sample, df in forStats:
                if sample == "A":
                    treat1 = df["Intensity"]
                    samples.append(treat1)
                elif sample == "B":
                    treat2 = df["Intensity"]
                    samples.append(treat2)
                elif sample == "C":
                    treat3 = df["Intensity"]
                    samples.append(treat3)
                elif sample == "D":
                    treat4 = df["Intensity"]
                    samples.append(treat4)
            f, p = stats.f_oneway(treat1, treat2, treat3, treat4)
            if p < 0.05:
                # test whether colimitation is different
                toControls, toControlp = stats.ttest_ind(treat1, treat4)
                toPs, toPp = stats.ttest_ind(treat2, treat4)
                toFes, toFep = stats.ttest_ind(treat3, treat4)
                if toPp and toFep < 0.05:
                    outCsv2.writerow([accession, sequence, f, p, toControlp, toPp, toFep])
                
                
                
        new_df_grouped = new_df.groupby("Sample")["Intensity"].agg([np.mean, np.std])

        
        try:
            intensityColimited =  new_df_grouped.loc["D"]['mean']
            for group in new_df_grouped.index:
                stdev = new_df_grouped.loc[group]['std']
                mean = new_df_grouped.loc[group]['mean']
                ratio = intensityColimited / mean
                # test the ratios
                # new_df_grouped['ratioTest'] = (ratio >=2 or ratio <= 0.5)
                if ratio >=2 or ratio <= 0.5:
                    outCsv.writerow([accession, description, sequence, group, mean, stdev, ratio])

        except KeyError:
            continue

def makePlots(directory, description):
    import os
    import csv
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    
    dataFile = os.path.join(directory, description + "significantDifferences.csv")
    df = pd.read_csv(dataFile)
    df_grouped = df.groupby("Sequence")
    sampleDict = {"A": "control", "B": "-P", "C": "-Fe", "D": "-P/-Fe"}
    
    psmsListFile = os.path.join(directory, description + "matchingpsms.csv")
    psms = pd.read_csv(psmsListFile, names = ["Accession", "Description", "Sequence", "Intensity", "Modifications", "Sample"])
    groupedBySequence = psms.groupby("Sample")["Intensity"].agg([np.mean, np.std])
    
    numplots = len(df_grouped)
    fig, ax = plt.subplots(numplots/2, 2, figsize = (10,numplots*2), dpi=20)
    a = 0
    b = 0
    
    for group, new_df in df_grouped:
        try:
            plot = groupedBySequence.loc[group].unstack()['mean'].plot(kind = "bar", ax = ax[a,b], yerr = groupedBySequence.loc[group, 'std'])
            plot.set_title(group)
            a += 1
        except IndexError:
            b += 1
            a = 0
                
    fig.tight_layout()
    outFig = os.path.join(directory, description +  ".png") # save to a new, numbered, figure
    plt.savefig(outFig)
    plt.close()
    

def main(directory, keywords, description):
    convertCSV(directory)
    getProteinsFromKeyword(directory, keywords, description)
    getPeptides(directory, description)
    makeComparisons (directory, description)
    #makePlots(directory, description)


if __name__ == "__main__":
    main('/Users/Noelle/Documents/tricho_colimitation/HD/', ["kinase"], "all_kinase")