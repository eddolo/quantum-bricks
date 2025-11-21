def apply_family_calibration(family, elems_set, Eg):
    if Eg <= 0:
        return Eg

    # keep these untouched
    if family in {"SI", "ZNO", "TIO2"}:
        return Eg
    if family == "GE":
        return Eg

    # elemental / simple
    if family == "DIAMOND":
        return 5.5

    # nitrides
    if family == "ALN":
        return 6.2
    if family == "GAN":
        return 3.4

    # oxides
    if family == "AL2O3":
        return 8.8
    if family == "MGO":
        return 7.8

    # II-VI already existing
    if family == "CDS":
        return 2.4
    if family == "CDSE":
        return 1.75

    # MoS2
    if family == "MOS2":
        return 1.2

    # oxide perovskites
    if family == "ABO3_PEROVSKITE":
        if Eg < 1.2:
            return 2.6
        elif Eg > 2.0:
            return 3.3
        else:
            return 3.0

    # halide perovskites
    if family == "ABX3_HALIDE":
        target = 1.6
        scale = target / Eg
        return max(0.4, min(2.0, Eg * scale))

    # -------------------------
    # NEW FAMILY CALIBRATIONS
    # -------------------------

    # III–V
    if family == "GAAS":
        return 1.42
    if family == "INP":
        return 1.35
    if family == "INAS":
        return 0.36
    if family == "ALAS":
        return 2.16
    if family == "GAP":
        return 2.3

    # II–VI extended
    if family == "ZNS":
        return 3.7
    if family == "ZNSE":
        return 2.7
    if family == "ZNTE":
        return 2.2
    if family == "CDTE":
        return 1.44

    # transition metal monoxides
    if family == "FEO":
        return 2.4
    if family == "COO":
        return 2.7
    if family == "NIO":
        return 3.7
    if family == "CUO":
        return 1.4
    if family == "MNO":
        return 2.0

    # SnX / PbX chalcogenides
    if family == "SNS":
        return 1.1
    if family == "SNSE":
        return 1.0
    if family == "SNTE":
        return 0.3
    if family == "PBS":
        return 0.4
    if family == "PBSE":
        return 0.27
    if family == "PBTE":
        return 0.31

    # simple halide salts
    if family == "NACL":
        return 8.5
    if family == "KCL":
        return 8.3
    if family == "CSCL":
        return 8.6
    if family == "LIF":
        return 14.2
    if family == "NAF":
        return 11.5

    return Eg
