#!/usr/bin/env python3

import os
import sys
import openpyxl


def Conv(sd, d=0):
    try:
        rt = int(sd)
    except Exception:
        rt = 0
    return str(rt)


def SeqChange(Ref, Alt):
    if (len(Ref) == 0) or (len(Alt) == 0):
        rt = ["", ""]
    elif len(Ref) < len(Alt):
        rt = [Ref[:1], Alt.rstrip(Ref[1:])]
    elif len(Ref) > len(Alt):
        rt = [Ref.rstrip(Alt[1:]), Alt[:1]]
    else:  # len(Ref) == len(Alt)
        rt = [Ref, Alt]
    return rt


class RAVclass:
    def __init__(self):
        # WHO DB for comparison
        script_dir = os.path.dirname(os.path.abspath(__file__))
        self.whofile = os.path.join(script_dir, "WHO_DB_V2_20241104.xlsx")
        
        self.varfile = ""
        self.wholist = []
        self.varlist = []

        # Indexes for fast lookup
        self.aa_index = {}   # (Gene, AminoAcid) -> WHO entry with highest Orderlev
        self.nuc_index = {}  # (StartPs, Ref, Alt) -> WHO entry

        self.get_whodata()

    def set_varfile(self, vfile):
        self.varfile = vfile

    def get_whodata(self):
        if not os.path.exists(self.whofile):
            print(self.whofile + " file does not exist")
            return False

        print("WHO DB Loading")
        # Identify required column index positions
        try:
            idx = {}
            cnt = 0

            xlsx = openpyxl.load_workbook(self.whofile, data_only=True)
            sheet = xlsx.worksheets[0]  # Sheet name: Annotation
            xlsx.close()

            for col in sheet["1"]:  # header row
                if col.value == "Target_Gene":
                    idx["Target_Gene"] = col.column - 1
                    cnt += 1
                elif col.value == "Accession_ID":
                    idx["Accession_ID"] = col.column - 1
                    cnt += 1
                elif col.value == "Start":
                    idx["Start"] = col.column - 1
                    cnt += 1
                elif col.value == "End":
                    idx["End"] = col.column - 1
                    cnt += 1
                elif col.value == "Ref":
                    idx["Ref"] = col.column - 1
                    cnt += 1
                elif col.value == "Alt":
                    idx["Alt"] = col.column - 1
                    cnt += 1
                elif col.value == "Variant_Type":
                    idx["Variant_Type"] = col.column - 1
                    cnt += 1
                elif col.value == "HGVS":
                    idx["HGVS"] = col.column - 1
                    cnt += 1
                elif col.value == "HGVS_AA":
                    idx["HGVS_AA"] = col.column - 1
                    cnt += 1
                elif col.value == "Tier":
                    idx["Tier"] = col.column - 1
                    cnt += 1
                elif col.value == "AMK":
                    idx["AMK"] = col.column - 1
                    cnt += 1
                elif col.value == "BDQ":
                    idx["BDQ"] = col.column - 1
                    cnt += 1
                elif col.value == "CAP":
                    idx["CAP"] = col.column - 1
                    cnt += 1
                elif col.value == "CFZ":
                    idx["CFZ"] = col.column - 1
                    cnt += 1
                elif col.value == "DLM":
                    idx["DLM"] = col.column - 1
                    cnt += 1
                elif col.value == "EMB":
                    idx["EMB"] = col.column - 1
                    cnt += 1
                elif col.value == "ETO":
                    idx["ETO"] = col.column - 1
                    cnt += 1
                elif col.value == "INH":
                    idx["INH"] = col.column - 1
                    cnt += 1
                elif col.value == "KAN":
                    idx["KAN"] = col.column - 1
                    cnt += 1
                elif col.value == "LFX":
                    idx["LFX"] = col.column - 1
                    cnt += 1
                elif col.value == "LZD":
                    idx["LZD"] = col.column - 1
                    cnt += 1
                elif col.value == "MFX":
                    idx["MFX"] = col.column - 1
                    cnt += 1
                elif col.value == "PZA":
                    idx["PZA"] = col.column - 1
                    cnt += 1
                elif col.value == "RIF":
                    idx["RIF"] = col.column - 1
                    cnt += 1
                elif col.value == "STM":
                    idx["STM"] = col.column - 1
                    cnt += 1
                elif col.value == "Note":
                    idx["Note"] = col.column - 1
                    cnt += 1
                elif col.value == "References":
                    idx["References"] = col.column - 1
                    cnt += 1
                elif col.value == "Order_level":
                    idx["Order_level"] = col.column - 1
                    cnt += 1

            if cnt != 28:
                print("WHO file header error")
                return False

        except Exception as e:
            print(e)
            return

        # Load WHO DB entries
        self.wholist = []
        self.aa_index = {}
        self.nuc_index = {}

        for r, rows in enumerate(sheet.rows):
            if r == 0:
                continue

            sid = rows[idx["Target_Gene"]].value
            if not sid:
                continue

            try:
                dt = {}
                dt["Gene"] = rows[idx["Target_Gene"]].value
                dt["GeneID"] = rows[idx["Accession_ID"]].value
                dt["StartPs"] = Conv(rows[idx["Start"]].value)
                dt["EndPs"] = Conv(rows[idx["End"]].value)
                dt["Ref"] = rows[idx["Ref"]].value
                dt["Alt"] = rows[idx["Alt"]].value
                dt["VarType"] = rows[idx["Variant_Type"]].value
                dt["Necleotide"] = rows[idx["HGVS"]].value
                dt["AminoAcid"] = rows[idx["HGVS_AA"]].value
                dt["Tier"] = Conv(rows[idx["Tier"]].value)
                dt["AMK"] = Conv(rows[idx["AMK"]].value)
                dt["BDQ"] = Conv(rows[idx["BDQ"]].value)
                dt["CAP"] = Conv(rows[idx["CAP"]].value)
                dt["CFZ"] = Conv(rows[idx["CFZ"]].value)
                dt["DLM"] = Conv(rows[idx["DLM"]].value)
                dt["EMB"] = Conv(rows[idx["EMB"]].value)
                dt["ETO"] = Conv(rows[idx["ETO"]].value)
                dt["INH"] = Conv(rows[idx["INH"]].value)
                dt["KAN"] = Conv(rows[idx["KAN"]].value)
                dt["LFX"] = Conv(rows[idx["LFX"]].value)
                dt["LZD"] = Conv(rows[idx["LZD"]].value)
                dt["MFX"] = Conv(rows[idx["MFX"]].value)
                dt["PZA"] = Conv(rows[idx["PZA"]].value)
                dt["RIF"] = Conv(rows[idx["RIF"]].value)
                dt["STM"] = Conv(rows[idx["STM"]].value)
                dt["Reference"] = rows[idx["References"]].value

                # Orderlev used for priority; ensure numeric
                order_val = rows[idx["Order_level"]].value
                try:
                    dt["Orderlev"] = float(order_val) if order_val is not None else 0.0
                except Exception:
                    dt["Orderlev"] = 0.0

                # Keep only rows where at least one drug has annotation
                if (
                    dt["AMK"] != "0"
                    or dt["BDQ"] != "0"
                    or dt["CAP"] != "0"
                    or dt["CFZ"] != "0"
                    or dt["DLM"] != "0"
                    or dt["EMB"] != "0"
                    or dt["ETO"] != "0"
                    or dt["INH"] != "0"
                    or dt["KAN"] != "0"
                    or dt["LFX"] != "0"
                    or dt["LZD"] != "0"
                    or dt["MFX"] != "0"
                    or dt["PZA"] != "0"
                    or dt["RIF"] != "0"
                    or dt["STM"] != "0"
                ):
                    self.wholist.append(dt)

                    # Build AA-based index: choose highest Orderlev per (Gene, AminoAcid)
                    key_aa = (dt["Gene"], dt["AminoAcid"])
                    if key_aa[1] not in (None, "", ".", "-"):
                        prev = self.aa_index.get(key_aa)
                        if (prev is None) or (dt["Orderlev"] > prev["Orderlev"]):
                            self.aa_index[key_aa] = dt

                    # Build nucleotide-based index: (StartPs, Ref, Alt)
                    key_nuc = (dt["StartPs"], dt["Ref"], dt["Alt"])
                    self.nuc_index[key_nuc] = dt

            except Exception as e:
                print(e)

    def get_vardata(self):
        if not os.path.isfile(self.varfile):
            print(self.varfile + " file does not exist")
            return

        if len(self.wholist) == 0:
            return

        self.varlist = []

        with open(self.varfile, "rt") as f:
            for r, rows in enumerate(f):
                # optional simple progress print (every 1000 lines)
                if r % 1000 == 0:
                    print(f"data loading line {r}")

                newitem = rows.rstrip("\n").split("\t")
                if len(newitem) < 10:
                    continue

                newdp = newitem[9].strip().split(",")

                if (r > 0) and (len(newitem) == 10) and (len(newdp) == 4):
                    seq = SeqChange(newitem[1], newitem[2])  # match WHO DB format
                    vty = newitem[4].split(",")[0]
                    gen = newitem[5].split(",")[0]
                    gid = newitem[6].split(",")[0]
                    nuc = newitem[7].split(",")[0]
                    aac = newitem[8].split(",")[0]

                    try:
                        dep = (
                            int(newdp[0])
                            + int(newdp[1])
                            + int(newdp[2])
                            + int(newdp[3])
                        )
                        if dep > 0:
                            frq = ((int(newdp[2]) + int(newdp[3])) / dep) * 100
                        else:
                            frq = 0.0
                    except Exception:
                        dep = 0
                        frq = 0.0

                    mat = 0
                    w = [""] * 16

                    # 1. WHO DB Comparison (Gene, Amino Acid Change)
                    #    Overlap level = Amino acid (mat = 2)
                    if aac not in [".", "-", "", " "]:
                        key_aa = (gen, aac)
                        who_aa = self.aa_index.get(key_aa)
                        if who_aa is not None:
                            w[0] = who_aa["Tier"]
                            w[1] = who_aa["AMK"]
                            w[2] = who_aa["BDQ"]
                            w[3] = who_aa["CAP"]
                            w[4] = who_aa["CFZ"]
                            w[5] = who_aa["DLM"]
                            w[6] = who_aa["EMB"]
                            w[7] = who_aa["ETO"]
                            w[8] = who_aa["INH"]
                            w[9] = who_aa["KAN"]
                            w[10] = who_aa["LFX"]
                            w[11] = who_aa["LZD"]
                            w[12] = who_aa["MFX"]
                            w[13] = who_aa["PZA"]
                            w[14] = who_aa["RIF"]
                            w[15] = who_aa["STM"]
                            mat = 2

                    # 2. WHO DB Comparison (Position, Ref, Alt)
                    #    Overlap level = Nucleotide (mat = 1)
                    key_nuc = (str(newitem[0]), seq[0], seq[1])
                    who_nuc = self.nuc_index.get(key_nuc)
                    if who_nuc is not None:
                        w[0] = who_nuc["Tier"]
                        w[1] = who_nuc["AMK"]
                        w[2] = who_nuc["BDQ"]
                        w[3] = who_nuc["CAP"]
                        w[4] = who_nuc["CFZ"]
                        w[5] = who_nuc["DLM"]
                        w[6] = who_nuc["EMB"]
                        w[7] = who_nuc["ETO"]
                        w[8] = who_nuc["INH"]
                        w[9] = who_nuc["KAN"]
                        w[10] = who_nuc["LFX"]
                        w[11] = who_nuc["LZD"]
                        w[12] = who_nuc["MFX"]
                        w[13] = who_nuc["PZA"]
                        w[14] = who_nuc["RIF"]
                        w[15] = who_nuc["STM"]
                        mat = 1  # this overwrites amino-acid-only match, as in original code

                    if mat > 0:
                        var = {}
                        var["Position"] = newitem[0]
                        var["Ref"] = newitem[1]
                        var["Alt"] = newitem[2]
                        var["Depth"] = newitem[3]
                        var["VarAlleleFreq"] = str(frq)
                        var["VarType"] = vty
                        var["Gene"] = gen
                        var["GeneID"] = gid
                        var["Necleotide"] = nuc
                        var["AminoAcid"] = aac
                        var["Tier"] = w[0]
                        var["AMK"] = w[1]
                        var["BDQ"] = w[2]
                        var["CAP"] = w[3]
                        var["CFZ"] = w[4]
                        var["DLM"] = w[5]
                        var["EMB"] = w[6]
                        var["ETO"] = w[7]
                        var["INH"] = w[8]
                        var["KAN"] = w[9]
                        var["LFX"] = w[10]
                        var["LZD"] = w[11]
                        var["MFX"] = w[12]
                        var["PZA"] = w[13]
                        var["RIF"] = w[14]
                        var["STM"] = w[15]
                        if mat == 1:
                            var["OverlapLv"] = "Nucleotide"
                        elif mat == 2:
                            var["OverlapLv"] = "Amino acid"
                        elif mat == 3:
                            var["OverlapLv"] = "Both"
                        else:
                            var["OverlapLv"] = ""

                        self.varlist.append(var)

    def save(self):
        if len(self.varlist) == 0:
            return

        out_path = self.varfile + ".txt"
        if os.path.isfile(out_path):
            os.remove(out_path)

        with open(out_path, "wt") as f:
            f.write(
                "Position\tREF\tALT\tDepth\tVariant allele frequency\tVariant type\t"
                "Gene\tGene ID\tNucleotide change\tAmino Acid change\tAMK\tBDQ\tCAP\tCFZ\t"
                "DLM\tEMB\tETO\tINH\tKAN\tLFX\tLZD\tMFX\tPZA\tRIF\tSTM\tOverlap Level\n"
            )
            for rows in self.varlist:
                s = (
                    rows["Position"]
                    + "\t"
                    + rows["Ref"]
                    + "\t"
                    + rows["Alt"]
                    + "\t"
                    + rows["Depth"]
                    + "\t"
                    + rows["VarAlleleFreq"]
                    + "\t"
                    + rows["VarType"]
                    + "\t"
                    + rows["Gene"]
                    + "\t"
                    + rows["GeneID"]
                    + "\t"
                    + rows["Necleotide"]
                    + "\t"
                    + rows["AminoAcid"]
                    + "\t"
                    + rows["AMK"]
                    + "\t"
                    + rows["BDQ"]
                    + "\t"
                    + rows["CAP"]
                    + "\t"
                    + rows["CFZ"]
                    + "\t"
                    + rows["DLM"]
                    + "\t"
                    + rows["EMB"]
                    + "\t"
                    + rows["ETO"]
                    + "\t"
                    + rows["INH"]
                    + "\t"
                    + rows["KAN"]
                    + "\t"
                    + rows["LFX"]
                    + "\t"
                    + rows["LZD"]
                    + "\t"
                    + rows["MFX"]
                    + "\t"
                    + rows["PZA"]
                    + "\t"
                    + rows["RIF"]
                    + "\t"
                    + rows["STM"]
                    + "\t"
                    + rows["OverlapLv"]
                    + "\n"
                )
                f.write(s)


if __name__ == "__main__":
    print("====================================")
    print(" Comparison with WHO DB 2nd edition ")
    print("====================================")

    if len(sys.argv) != 2:
        print("Usage: python RAV_v1.1.py <directory_path>")
    elif not os.path.isdir(sys.argv[1]):
        print("Directory path does not exist.")
    else:
        directory_path = sys.argv[1]

        # Load WHO DB once and reuse for all TSV files
        rav = RAVclass()

        for filename in os.listdir(directory_path):
            if filename.endswith(".tsv"):
                full_path = os.path.join(directory_path, filename)
                print(f"Processing: {full_path}")
                rav.set_varfile(full_path)
                rav.get_vardata()
                rav.save()
                print(f"Processed: {filename}")

        print("Finish!")