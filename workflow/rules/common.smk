def get_h1_tag(wc):
    """Get the H1 tag for a sample."""
    return MANIFEST[wc.sm]["h1_tag"]


def get_h2_tag(wc):
    """Get the H2 tag for a sample."""
    return MANIFEST[wc.sm]["h2_tag"]


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