#
# this script loads in tree seqs outputted by slim
# exports individual locations and trait values to csv's
# then exports genotype arrays of causal mutations
# then wipes the causal mutations from the tree seq
# then sprinkles on neutral mutations
# then exports genotype arrays of neutral mutations
#

import tskit
import msprime
import numpy as np

def genotypeArrays(datapth):

    # load trees

    host_ts = tskit.load(datapth+'host-slim.trees')
    para_ts = tskit.load(datapth+'para-slim.trees')


    # save causal snp locations

    h_csl_snps = []
    for i in host_ts.sites():
        h_csl_snps.append(i.position)
    np.savetxt(datapth+'h-csl-snps.csv', h_csl_snps, delimiter=",")

    p_csl_snps = []
    for i in para_ts.sites():
        p_csl_snps.append(i.position)
    np.savetxt(datapth+'p-csl-snps.csv', p_csl_snps, delimiter=",")


    # export locations and trait values

    indivlist = []
    indivnames = []
    with open(datapth+'host.txt', "w") as indfile:
        indfile.writelines(",".join(["vcf_label"]
                                + ["x", "y", "z"]) + "\n") # (x,y)=coords, z=trait
        for i in host_ts.individuals():
            indivlist.append(i.id)
            vcf_label = f"tsk_{i.id}"
            indivnames.append(vcf_label)
            data = [vcf_label, str(i.location[0]), str(i.location[1]), str(i.location[2])]
            indfile.writelines(",".join(data) + "\n")

    indivlist = []
    indivnames = []
    with open(datapth+'para.txt', "w") as indfile:
        indfile.writelines(",".join(["vcf_label"]
                                + ["x", "y", "z"]) + "\n") # (x,y)=coords, z=trait
        for i in para_ts.individuals():
            indivlist.append(i.id)
            vcf_label = f"tsk_{i.id}"
            indivnames.append(vcf_label)
            data = [vcf_label, str(i.location[0]), str(i.location[1]), str(i.location[2])]
            indfile.writelines(",".join(data) + "\n")

    # export individual allele freqs and genotype matrices of causal mutations

    hgm = host_ts.genotype_matrix()
    hgm1 = hgm[:,0:host_ts.num_individuals]
    hgm2 = hgm[:,host_ts.num_individuals:host_ts.num_samples]
    hfrq = (hgm1+hgm2)/2
    np.save(datapth+'hfrq-causal',hfrq)
    hga = np.reshape(np.array([hgm1,hgm2]),(host_ts.num_sites,host_ts.num_individuals,2))
    np.save(datapth+'hga-causal',hga)

    pgm = para_ts.genotype_matrix()
    pgm1 = pgm[:,0:para_ts.num_individuals]
    pgm2 = pgm[:,para_ts.num_individuals:para_ts.num_samples]
    pfrq = (pgm1+pgm2)/2
    np.save(datapth+'pfrq-causal',pfrq)
    pga = np.reshape(np.array([pgm1,pgm2]),(para_ts.num_sites,para_ts.num_individuals,2))
    np.save(datapth+'pga-causal',pga)

    # sprinkle on neutral mutations

    mu = 1e-11 # neutral mutations occur an order of magnitude faster than causal

    htables = host_ts.tables
    htables.mutations.clear()
    htables.sites.clear()
    hts = htables.tree_sequence()

    ptables = para_ts.tables
    ptables.mutations.clear()
    ptables.sites.clear()
    pts = ptables.tree_sequence()

    hmut = msprime.SLiMMutationModel(type=1)
    hts = msprime.sim_mutations(hts,rate=mu,model=hmut,keep=True)

    pmut = msprime.SLiMMutationModel(type=1)
    pts = msprime.sim_mutations(pts,rate=mu,model=pmut,keep=True)

    # save neutral snp locations

    h_ntl_snps = []
    for i in hts.sites():
        h_ntl_snps.append(i.position)
    np.savetxt(datapth+'h-ntl-snps.csv', h_ntl_snps, delimiter=",")

    p_ntl_snps = []
    for i in pts.sites():
        p_ntl_snps.append(i.position)
    np.savetxt(datapth+'p-ntl-snps.csv', p_ntl_snps, delimiter=",")

    # export individual allele freqs and genotype matrices of neutral mutations

    hgm = hts.genotype_matrix()
    hgm1 = hgm[:,0:hts.num_individuals]
    hgm2 = hgm[:,hts.num_individuals:hts.num_samples]
    hfrq = (hgm1+hgm2)/2
    np.save(datapth+'hfrq-neutrl',hfrq)
    hga = np.reshape(np.array([hgm1,hgm2]),(hts.num_sites,hts.num_individuals,2))
    np.save(datapth+'hga-neutrl',hga)

    pgm = pts.genotype_matrix()
    pgm1 = pgm[:,0:pts.num_individuals]
    pgm2 = pgm[:,pts.num_individuals:pts.num_samples]
    pfrq = (pgm1+pgm2)/2
    np.save(datapth+'pfrq-neutrl',pfrq)
    pga = np.reshape(np.array([pgm1,pgm2]),(pts.num_sites,pts.num_individuals,2))
    np.save(datapth+'pga-neutrl',pga)