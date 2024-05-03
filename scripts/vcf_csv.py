import sys
import vcf

vcf_file = sys.argv[1] #pass vcf file
output_csv = sys.argv[2] #file to write

vcf_reader = vcf.Reader(open(vcf_file, 'r')) #open VCF file
sample_names = vcf_reader.samples #extract sample names

# Open CSV file for writing
with open(output_csv, 'w') as csv_file:
    csv_file.write("CHROM;POS;ID;REF;ALT;QUAL;FILTER;FORMAT;|Annotation|Annotation_Impact|Gene_Name|Gene_ID|Feature_Type|;" + ";".join(sample_names) + "\n") #write header

    for record in vcf_reader: #write data
        info=record.INFO['ANN']
        info_string = "".join(str(element) for element in info)
        var_impact = info_string.split('|')[1:6]
        var_impact_string = "|".join(str(element) for element in var_impact)
        csv_file.write(f"{record.CHROM};{record.POS};{record.ID};{record.REF};{record.ALT[0] if record.ALT else '.'};{record.QUAL or '.'};{record.FILTER[0] if record.FILTER else '.'};DP:RO:AO:QA;{var_impact_string};")
        for sample in record.samples: #write genotype information for each sample
            if sample.data.GT:

                csv_file.write(f"{sample.data.DP}:{sample.data.RO}:{sample.data.AO}:{sample.data.QA};")

        csv_file.write("\n")
print("Conversion completed.")