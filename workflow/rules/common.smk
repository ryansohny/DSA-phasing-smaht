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
