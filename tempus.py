import vcf
import time
import urllib.request

def covered(variant):
    #be nice to the server :)
    time.sleep(0.3)
    with urllib.request.urlopen('http://exac.hms.harvard.edu/rest/variant/any_covered/'+variant) as response:
        buff = response.read()
        return('true' in str(buff))

def fetch_allele_freq(variant):
    # be nice to the server :)
    time.sleep(0.3)
    with urllib.request.urlopen('http://exac.hms.harvard.edu/rest/variant/variant/'+variant) as response:
        buff = response.read()
        dic=eval(buff)
        if('allele_freq' in dic):
            return(dic['allele_freq'])
        else:
            return(0)

def fetch_consequences(variant):
    # be nice to the server :)
    time.sleep(0.3)
    with urllib.request.urlopen('http://exac.hms.harvard.edu/rest/variant/consequences/'+variant) as response:
        buff = eval(response.read())
        return(buff)


out=open('tempus.txt','w')
vcf_reader = vcf.Reader(filename='t.vcf')
print("VariantID\tDepth\tSupportingReads\tSupportingReadsPercentage\tExAC_VariantFrequency\tExAC_Consequences")
out.write("VariantID\tDepth\tSupportingReads\tSupportingReadsPercentage\tExAC_VariantFrequency\tExAC_Consequences\n")
count=0;
for record in vcf_reader:
    allels=len(record.ALT)
    #debugging insertion
    if(count>100):
        break
    count+=1
    for i in range(allels):
        variantID = str(record.CHROM)+'-'+str(record.POS)+'-'+str(record.REF)+'-'+str(record.ALT[i])
        percent = (record.INFO['AO'][i]*100)/record.INFO['DP']
        known = covered(variantID)
        if(known== True):
            freq = fetch_allele_freq(variantID)
            if(freq>0):
                cons = fetch_consequences(variantID)
                print("%s\t%d\t%d\t%d\t%f\t%s" % (variantID, record.INFO['DP'], record.INFO['AO'][i], percent, freq, list(cons)))
                out.write("%s\t%d\t%d\t%d\t%f\t%s\n" % (variantID, record.INFO['DP'], record.INFO['AO'][i], percent, freq, list(cons)))
            else:
                print("%s\t%d\t%d\t%d\tincorrect call\tunclear" % (variantID, record.INFO['DP'], record.INFO['AO'][i], percent))
                out.write("%s\t%d\t%d\t%d\tincorrect call\tunclear\n" % (variantID, record.INFO['DP'], record.INFO['AO'][i], percent))
        else:
            print("%s\t%d\t%d\t%d\tunknown variant\tunknown consequences" % (variantID, record.INFO['DP'], record.INFO['AO'][i], percent))
            out.write("%s\t%d\t%d\t%d\tunknown variant\tunknown consequences\n" % (variantID, record.INFO['DP'], record.INFO['AO'][i], percent))