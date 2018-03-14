#!/usr/bin python3

""" Search using Accession number """



def SearchAcc(acc):

    """ Use genbank accession number to retrieve other identifiers for
    specified gene

    Input:  acc                                --- Genbank accession number 
    
    Output: (acc, genid, product, location)    --- Accession number, Gene id,
                                                   protein product, chromosomal location

    09.03.18    Original    By:JJS
    14.03.18                By:JJS
            
    """

       
    # checks if accession number is in dictionary
    # returns gene id, protein product and chrom location

    if acc in acc_dict:
        data={'acc':acc, 'genid':genid, 'product':product, 'location':location}
        return(data)
    else:
        pass


dummy = 'AB061209'
genid = 'MRPS28'
product = 'mitochondrial ribosomal protein s28'
location = '8q21.1-q21.2'

acc_dict = {'AB061209', 'AB098765' 'AC012345'}

with open("acc_output.txt", "w") as text_file:
    xmlTemplate = """<gene>
    <geneidentifiers>
        <accession>%(acc)s</accession>
        <genid>%(genid)s</genid>
        <product>%(product)s</product>
        <location>%(location)s</location>
    </geneidentifiers>
</gene>"""
    
    test = SearchAcc(dummy)
    print(xmlTemplate%test, file=text_file)

