import hail as hl
import os

hl.init(
    default_reference='GRCh38',
    tmp_dir='/data/tmp/',
    local_tmpdir='/data/tmp/'
    )

vcf_file = "/data7/ipsc/16p12_1_del/str_calls/data/filtered_data/CR001/CR001.vcf"
matrix_table = '/data7/ipsc/16p12_1_del/str_calls/data/matrix_tables/CR001.mt'


def create_hail_matrix(vcf_file, matrix_table):
    os.makedirs(os.path.dirname(matrix_table), exist_ok=True)
    hl.import_vcf(vcf_file, array_elements_required=False).write(matrix_table, overwrite=True)
    return

if __name__ == "__main__":
    vcf_file = "/data7/ipsc/16p12_1_del/str_calls/data/filtered_data/CR001/CR001.vcf"
    matrix_table = '/data7/ipsc/16p12_1_del/str_calls/data/matrix_tables/CR001.mt'
    create_hail_matrix(vcf_file, matrix_table)
