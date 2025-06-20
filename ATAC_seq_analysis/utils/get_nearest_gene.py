import os
import pybedtools

#### GLOBALS ####
CURRENT_DIR_PATH = os.path.dirname(os.path.abspath(__file__))
pybedtools.helpers.set_tempdir(CURRENT_DIR_PATH)


def get_closest_gene(master_file, gtf_file, save_file):
    a = pybedtools.BedTool(master_file)
    # get the closest feature in 'other.bed' on the same strand
    gtf_bed = pybedtools.BedTool(gtf_file).sort()
    b = a.closest(gtf_bed, d=True)
    b.moveto(save_file)
    return


if __name__ == "__main__":
    master_file = "npc_atac_peaks.bed"
    gtf_file = "gencode.v44.basic.annotation.gtf"
    save_file = "npc_closest.bed"

    get_closest_gene(master_file, gtf_file, save_file)
    pybedtools.helpers.cleanup(verbose=False, remove_all=False)
