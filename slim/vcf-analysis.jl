using GeneticVariation
using VCFTools

reader = VCF.Reader(open("/home/bb/gits/genomic-sign-coev-cont-sp/SLiM/host.vcf", "r"))
for record in reader
    # do something
end
close(reader)

Base.find(header(reader), "FORMAT")


fh = openvcf("/home/bb/gits/genomic-sign-coev-cont-sp/SLiM/host.vcf", "r");
for l in 1:4
    println(readline(fh))
end
close(fh)



nrecords("/home/bb/gits/genomic-sign-coev-cont-sp/SLiM/host.vcf")

nsamples("/home/bb/gits/genomic-sign-coev-cont-sp/SLiM/host.vcf")

@time A = convert_gt(Int8, "/home/bb/gits/genomic-sign-coev-cont-sp/SLiM/host.vcf"; model = :additive, impute = false, center = false, scale = false)

@time records, samples, lines, missings_by_sample, missings_by_record, 
    maf_by_record, minorallele_by_record = gtstats("/home/bb/gits/genomic-sign-coev-cont-sp/SLiM/host.vcf");