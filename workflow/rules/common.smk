def is_ont(wc):
    """True iff the sample's platform is ONT."""
    return MANIFEST[wc.sm]["platform"] == "ont"


def is_illumina(wc):
    """True iff the sample's platform is Illumina."""
    return MANIFEST[wc.sm]["platform"] == "illumina"


def is_fiberseq(wc):
    """True iff the sample is Fiber-seq (Hia5-treated)."""
    return MANIFEST[wc.sm]["fiber-seq"]


def get_h1_tag(wc):
    """Get the H1 tag for a sample."""
    h1 = f"'{MANIFEST[wc.sm]['h1_tag']}'"
    return h1


def get_h2_tag(wc):
    """Get the H2 tag for a sample."""
    h2 = f"'{MANIFEST[wc.sm]['h2_tag']}'"
    return h2


def get_dsa(wc):
    """Get the DSA for a sample."""
    return MANIFEST[wc.sm]["dsa"]


def get_mm2_preset(wc):
    """minimap2 -ax preset for the initial DSA alignment, chosen by platform.

    PacBio uses 'lr:hqae', ONT uses 'lr:hq'. An explicit `mm2_preset` in the
    config overrides the per-platform default for every sample. (Shared-reference
    realignment is a separate process and always uses 'lr:hq'.)
    """
    if config.get("mm2_preset"):
        return config["mm2_preset"]
    return {"pacbio": "'lr:hqae'", "ont": "'lr:hq'"}[MANIFEST[wc.sm]["platform"]]


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


def get_pre_haplotag_input(wc):
    """BAM source for haplotag_and_sort: minimap2 BAM (long-read) or bwa BAM (Illumina)."""
    if is_illumina(wc):
        return f"temp/{wc.sm}.{wc.file_idx}.illumina.bam"
    return f"temp/{wc.sm}.{wc.file_idx}.bam"


def get_fire_input(wc):
    """Get the alignment file input for FIRE - either modkit BAM or haplotag_and_sort CRAM."""
    if is_ont(wc):
        # ONT: FIRE runs after modkit (BAM output)
        return f"temp/{wc.sm}.{wc.file_idx}.modkit.dsa.bam"
    else:
        # PacBio: FIRE runs after haplotag_and_sort (CRAM output)
        return f"temp/{wc.sm}.{wc.file_idx}.dsa.cram"


def get_crams_to_merge(wc):
    """Pick the per-file CRAM source for merge_sample based on Fiber-seq status."""
    if is_fiberseq(wc):
        template = "temp/{{sm}}.{file_idx}.fire.dsa.cram"
    else:
        template = "temp/{{sm}}.{file_idx}.dsa.cram"
    return expand(template, file_idx=get_file_indices(wc.sm))


def get_crais_to_merge(wc):
    """Get all CRAI index files for a sample."""
    return [f"{cram}.crai" for cram in get_crams_to_merge(wc)]


def get_final_cram(wc):
    """Get the final CRAM file for a sample (merged output for both PacBio and ONT)."""
    return f"results/{wc.sm}.dsa.cram"


def bam_header_sm_settings(wc):
    if config.get("set-sm", False):
        return f" --sample {wc.sm} "
    else:
        return ""


def read_assignment_result():
    if config.get("keep_read_assignments", False):
        return "results/{sm}/{sm}.{file_idx}.assignments.tsv.gz"
    return temp("temp/{sm}.{file_idx}.assignments.tsv.gz")


def samples_by(predicate):
    """Return a regex alternation of sample IDs for which predicate(sm) is True."""
    matching = [sm for sm in SMs if predicate(sm)]
    if not matching:
        # Use a regex that matches nothing so the rule never fires.
        return r"(?!.*)"
    return "|".join(matching)


FIBERSEQ_SMS = samples_by(lambda sm: MANIFEST[sm]["fiber-seq"])
LONGREAD_NOFS_SMS = samples_by(
    lambda sm: not MANIFEST[sm]["fiber-seq"] and MANIFEST[sm]["platform"] != "illumina"
)
ILLUMINA_SMS = samples_by(lambda sm: MANIFEST[sm]["platform"] == "illumina")
