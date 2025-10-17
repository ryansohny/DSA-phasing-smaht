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
    """Get the BAM file(s) for a sample."""
    bam_files = MANIFEST[wc.sm]["bam"]
    # For rules that need a single BAM, return the first one
    # or the specific one based on file_idx
    if hasattr(wc, "file_idx"):
        return bam_files[int(wc.file_idx)]
    return bam_files


def get_all_bams(wc):
    """Get all BAM files for a sample."""
    return MANIFEST[wc.sm]["bam"]


def get_num_files(wc):
    """Get the number of input files for a sample."""
    return len(MANIFEST[wc.sm]["bam"])


def get_file_indices(sm):
    """Get list of file indices for a sample."""
    return list(range(len(MANIFEST[sm]["bam"])))


def get_crams_to_merge(wc):
    """Get all CRAM files for a sample - either haplotagged or modkit based on ONT tag."""
    is_ont = MANIFEST[wc.sm].get("ont", False) or config.get("ont", False)
    if is_ont:
        # ONT: merge after modkit
        return expand(
            "temp/{{sm}}.{file_idx}.modkit.dsa.cram", file_idx=get_file_indices(wc.sm)
        )
    else:
        # PacBio: merge after haplotag_and_sort
        return expand(
            "temp/{{sm}}.{file_idx}.dsa.cram", file_idx=get_file_indices(wc.sm)
        )


def get_final_cram(wc):
    """Get the final CRAM file for a sample (merged output for both PacBio and ONT)."""
    return f"results/{wc.sm}.dsa.cram"


def bam_header_sm_settings(wc):
    if config.get("set-sm", False):
        return f" --sample {wc.sm} "
    else:
        return ""
