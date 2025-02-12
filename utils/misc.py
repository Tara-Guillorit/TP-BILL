
def parse_range(range_str):
    """Convert a range string (e.g. '1-5') to a list of integers."""
    # Split the string by the dash '-'
    start, end = range_str.split('-')
    start, end = int(start), int(end)
    
    # Generate a list of numbers from start to end (inclusive)
    return list(range(start, end + 1))

def build_vcf_path(sample, it):
    root = "/students/BILL/2025-BILL/"
    path_iter = {
        15: "P15/KHV-U_trunc/",
        30: "P30/KHV-U_trunc/",
        50: "P50/KHV-U_trunc/",
        65: "P65_new/vcf", 
        90: "P90/vcf/"
    }

    if it == 65:
        return root + path_iter[it] + f"P{it}.{sample}.trimed1000.sv_sniffles.vcf"
    
    return root + path_iter[it] + f"P{it}-{sample}.trimed1000.sv_sniffles.vcf"\

def build_vcf_path_test(sample, it):
    root = "data/"
    path_iter = {
        15: "data-p15/",
        30: "data-p30/",
        50: "data-p50/",
        65: "data-p65/", 
        90: "data-p90/"
    }
    
    return root + path_iter[it] + f"P{it}-{sample}.trimed1000.sv_sniffles.vcf"