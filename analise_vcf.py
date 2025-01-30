#il faut intaller la librairie suivante
#pip install pyvcf3
import vcf



vcf_reader = vcf.Reader(open('test.vcf','r'))
for record in vcf_reader:
    print (record)
