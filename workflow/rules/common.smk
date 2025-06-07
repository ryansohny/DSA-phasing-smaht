def get_h1_tag(wc):
    """Get the H1 tag for a sample."""
    h1 = f"'{MANIFEST[wc.sm]["h1_tag"]}'"
    return h1


def get_h2_tag(wc):
    """Get the H2 tag for a sample."""
    h2 = f"'{MANIFEST[wc.sm]["h2_tag"]}'"
    return h2


def get_dsa(wc):
    """Get the DSA for a sample."""
    return MANIFEST[wc.sm]["dsa"]


def get_bam(wc):
    """Get the BAM file for a sample."""
    return MANIFEST[wc.sm]["bam"]


def get_final_cram(wc):
    """Get the final CRAM file for a sample."""
    if config.get("ont", False):
        return f"results/{wc.sm}.modkit.dsa.cram"
    else:
        return f"results/{wc.sm}.dsa.cram"

def bam_header_sm_settings(wc):
    if config.get("set-sm", False):
        return f" --sample {wc.sm} "
    else:
        return ""