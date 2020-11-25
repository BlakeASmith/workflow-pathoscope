def find_sam_align_score(fields):
    """
    Find the Bowtie2 alignment score for the given split line (``fields``).

    Searches the SAM fields for the ``AS:i`` substring and extracts the Bowtie2-specific alignment score. This will not
    work for other aligners.

    :param fields: a line that has been split on "\t"
    :type fields: list

    :return: the alignment score
    :rtype: float

    """
    read_length = float(len(fields[9]))

    for field in fields:
        if field.startswith("AS:i:"):
            a_score = int(field[5:])
            return a_score + read_length

    raise ValueError("Could not find alignment score")
