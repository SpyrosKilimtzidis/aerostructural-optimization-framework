import numpy as np
import itertools

def select_laminates(xn, hply, compid):
    def fullfact(levels):
        return np.array(list(itertools.product(*[range(1, l + 1) for l in levels])))

    # Baseline layup generation
    dFF1 = fullfact([20, 20, 20])
    percent1 = []
    tlpud1 = []

    for row in dFF1:
        blinelp1 = ([45] * row[2] + [0] * row[0] + [-45] * row[2] + [90] * row[2])
        blinelp = blinelp1 + blinelp1[::-1]

        perc0 = blinelp.count(0)
        perc90 = blinelp.count(90)
        perc45 = blinelp.count(45)
        percm45 = blinelp.count(-45)
        total = len(blinelp)

        p0 = perc0 / total * 100
        p90 = perc90 / total * 100
        p45m45 = (perc45 + percm45) / total * 100
        valid = 1 if min(p0, p90, p45m45) >= 10 else 0
        percent1.append([p0, p90, p45m45, valid])

        tlpud1.append([hply[0] * total, hply[1] * total, 1])

    percent1 = np.array(percent1)
    tlpud1 = np.array(tlpud1)
    udf1 = np.hstack((dFF1, tlpud1, percent1))
    udfinal2 = udf1[np.argsort(udf1[:, 3])]  # sort by thickness

    laminate = []
    fthick = []

    for jk in range(len(xn)):
        # Select columns based on component ID
        if compid[jk] == 1:
            udfinal1 = np.hstack((udfinal2[:, 0:4], udfinal2[:, 5:]))
        else:
            udfinal1 = np.hstack((udfinal2[:, 0:3], udfinal2[:, 4:]))

        # Filter valid laminates
        mask_valid = (udfinal1[:, 3] >= xn[jk]) & (udfinal1[:, 8] == 1)
        udthick = udfinal1[mask_valid]

        # Desired vector for Euclidean distance
        desper = np.linalg.norm([62.5, 12.5, 25, xn[jk]]) if compid[jk] == 1 else np.linalg.norm([25, 25, 50, xn[jk]])
        diffs = []
        for row in udthick:
            diff = np.linalg.norm([row[5], row[6], row[7], row[3]]) - desper
            diffs.append(abs(diff))
        udthick = np.hstack((udthick, np.array(diffs).reshape(-1, 1)))

        # Target percentage filters
        if compid[jk] == 1:
            mask_target = (udthick[:, 5] >= 45) & (udthick[:, 6] >= 10) & (udthick[:, 7] >= 10)
        else:
            mask_target = (udthick[:, 5] >= 20) & (udthick[:, 6] >= 20) & (udthick[:, 7] >= 45)
        udthick2 = udthick[mask_target]

        # Thickness tolerance
        udthick3 = udthick2[udthick2[:, 3] - xn[jk] <= 2e-3]

        # Select closest laminate by thickness
        finalam1 = udthick3[np.argsort(udthick3[:, 3])]
        xtr = finalam1[0, 0:3].astype(int)

        # Construct laminate layup
        layup = ([45] * xtr[2] + [0] * xtr[0] + [-45] * xtr[2] + [90] * xtr[1] +
                 [90] * xtr[1] + [-45] * xtr[2] + [0] * xtr[0] + [45] * xtr[2])
        laminate.append(layup)
        fthick.append(finalam1[0, 3])

    return laminate, np.array(fthick)

def select_laminates_single(xn, hply, compid):
    laminate = []
    fthick = []
    fthick1 = []
    for jk in range(len(xn)):
        if compid[jk]==1:
            hply1 = hply[0]
        else:
            hply1 = hply[1]
        if xn[jk] <= hply1:
            layup = [0]
            fthick1 = hply1
        elif xn[jk] <= 2 * hply1:
            layup = [0, 90]
            fthick1 = 2*hply1
        else:
            layup = [0, 90, 0]
            fthick1 = 3*hply1
        laminate.append(layup)
        fthick.append(fthick1)
    return laminate, fthick