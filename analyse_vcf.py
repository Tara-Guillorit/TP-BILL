#il faut intaller la librairie suivante
#pip install pyvcf3
#https://github.com/Tara-Guillorit/TP-BILL
import vcf



vcf_reader = vcf.Reader(open('test.vcf','r'))
for record in vcf_reader:
    print (record)
