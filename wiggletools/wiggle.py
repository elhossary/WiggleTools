import pandas as pd
import os
import csv
import numpy as np
import re


class Wiggle:
    def __init__(self, file_path, chrom_sizes):
        self.file_path = file_path
        self.chrom_sizes = chrom_sizes
        self.wiggle_df = None
        self.raw_data = []
        self.parse()
        self.orientation = None

    def parse(self):
        current_wiggle_meta = {}
        with open(self.file_path, "r") as raw_file:
            print(f"==> Loading file: {os.path.basename(self.file_path)}")
            file_header, all_contents = self._parse_wiggle_str(raw_file.read())
            current_wiggle_meta = self.parse_wiggle_header(file_header, current_wiggle_meta)
            tmp_dict = {}
            for content_header, content in all_contents.items():
                if "-" in content:
                    self.orientation = "r"
                else:
                    self.orientation = "f"
                for i in content.split("\n"):
                    line_split = i.split(" ")
                    if len(line_split) == 2:
                        tmp_dict[int(line_split[0])] = float(line_split[1])
                current_wiggle_meta = self.parse_wiggle_header(content_header, current_wiggle_meta)
                self.raw_data.append({"track_type": current_wiggle_meta["track_type"],
                                      "track_name": current_wiggle_meta["track_name"],
                                      "variableStep_chrom": current_wiggle_meta["variableStep_chrom"],
                                      "variableStep_span": current_wiggle_meta["variableStep_span"],
                                      "data": tmp_dict.copy()})


    @staticmethod
    def _parse_wiggle_str(in_str):
        ret_dict = {}
        header_text = in_str.split("\n", maxsplit=1)[0]
        in_str = in_str.replace(header_text + "\n", "")
        all_headers = re.findall(r'^.*chrom=.*$', in_str, flags=re.MULTILINE | re.IGNORECASE)
        splitters = ""
        for header in all_headers:
            splitters += header + "|"
        splitters = f"({splitters[:-1]})"
        split_str_list = re.split(rf"{splitters}", in_str, flags=re.MULTILINE | re.IGNORECASE)
        content_list = [i for i in split_str_list if i != '']
        for i in range(0, len(content_list), 2):
            ret_dict[content_list[i]] = content_list[i + 1]
        return header_text, ret_dict

    @staticmethod
    def parse_wiggle_header(line, current_wiggle_meta):
        if "type=" in line:
            current_wiggle_meta["track_type"] = line.split('type=')[-1].split(' ')[0].replace('\"', '')
        if "name=" in line:
            current_wiggle_meta["track_name"] = line.split('name=')[-1].replace('\n', '').replace('\"', '')
        if "chrom=" in line:
            current_wiggle_meta["variableStep_chrom"] = line.split('chrom=')[-1].split(' ')[0].replace('\"', '')
        if "span=" in line:
            current_wiggle_meta["variableStep_span"] = line.split('span=')[-1].replace('\n', '').replace('\"', '')
        return current_wiggle_meta

    def _to_dataframe(self, is_full=True):
        ignored_seqid = []
        wiggle_seqid = list(set([item["variableStep_chrom"] for item in self.raw_data]))
        # drop coverage for sequences not in chrom sizes
        self.raw_data = [item for item in self.raw_data
                         if item["variableStep_chrom"] in [x["seqid"] for x in self.chrom_sizes]]
        intersect_seqid = [item["variableStep_chrom"] for item in self.raw_data]
        ignored_seqid.extend([x["seqid"] + " (found in fasta not in wiggle)"
                              for x in self.chrom_sizes if x["seqid"] not in intersect_seqid])
        ignored_seqid.extend([x + " (found in wiggle not in fasta)" for x in wiggle_seqid if x not in intersect_seqid])
        if is_full:
            self._extend_raw_data_to_full_length_dataframe()
        else:
            self._extend_raw_data_to_min_length_dataframe()
        condition_name = self.wiggle_df["track_name"].unique().tolist()

        # Logging
        s = '\n     └── '
        if len(ignored_seqid) == 0:
            print(f"===> Parsed condition: {', '.join(condition_name)}\n"
                  f"     + Included sequence IDs:\n"
                  f"     └── {s.join(intersect_seqid)}")
        else:
            print(f"===> Parsed condition: {', '.join(condition_name)}\n"
                  f"     + Included sequence IDs:\n"
                  f"     └── {s.join(intersect_seqid)}\n"
                  f"     + Ignored sequence IDs:\n"
                  f"     └── {s.join(ignored_seqid)}")

    def _extend_raw_data_to_full_length_dataframe(self):
        extended_data = []
        for item in self.raw_data:
            meta_list = [item["track_type"],
                         item["track_name"],
                         item["variableStep_chrom"],
                         item["variableStep_span"]]
            # make backbone
            chrom_size = 0
            for chrom in self.chrom_sizes:
                if chrom["seqid"] == item["variableStep_chrom"]:
                    chrom_size = chrom["size"]
            tmp_list = [x for x in range(1, chrom_size + 1, 1)]
            # map the scores to backbone
            for index, val in enumerate(tmp_list):
                if val in item["data"].keys():
                    tmp_list[index] = meta_list + [val, item["data"][val]]
                else:
                    tmp_list[index] = meta_list + [val, 0.0]
            extended_data.extend(tmp_list)

        column_names = ["track_type", "track_name", "variableStep_chrom", "variableStep_span", "location", "score"]
        self.wiggle_df = pd.DataFrame(extended_data, columns=column_names)

    def _extend_raw_data_to_min_length_dataframe(self):
        extended_data = []
        for item in self.raw_data:
            meta_list = [item["track_type"],
                         item["track_name"],
                         item["variableStep_chrom"],
                         item["variableStep_span"]]
            for k, v in item["data"].items():
                extended_data.append(meta_list + [k, v])
        column_names = ["track_type", "track_name", "variableStep_chrom", "variableStep_span", "location", "score"]
        self.wiggle_df = pd.DataFrame(extended_data, columns=column_names)

    def get_wiggle(self, is_full=True):
        self._to_dataframe(is_full)
        return self.wiggle_df

    def write_wiggle(self, out_path, alt_wiggle_df=None):

        print(f"==> Writing file: {os.path.basename(out_path)}")
        if alt_wiggle_df is None:
            out_df = self.wiggle_df.copy()
        else:
            out_df = alt_wiggle_df.copy()
        # Shrink the dataframe
        out_df = out_df[out_df["score"] != 0]

        out_df["loc_score"] = out_df["location"].astype(str) + " " + out_df["score"].astype(str)
        out_df.drop(["location", "score"], axis=1, inplace=True)
        out_df = out_df.groupby(["track_type", "track_name", "variableStep_chrom", "variableStep_span"])["loc_score"]\
            .apply('\n'.join).reset_index()
        out_df["header"] = 'track type=' + out_df["track_type"] + ' name="' + out_df["track_name"] + '"'
        out_df.drop(["track_type", "track_name"], axis=1, inplace=True)
        out_df["sub_header_with_data"] = 'variableStep chrom=' + out_df["variableStep_chrom"] + \
                                         ' span=' + out_df["variableStep_span"] + '\n' + out_df["loc_score"]
        out_df.drop(["variableStep_chrom", "variableStep_span", "loc_score"], axis=1, inplace=True)
        out_df = out_df.groupby(["header"])["sub_header_with_data"].apply('\n'.join).reset_index()
        with open(out_path, 'w') as f:
            f.write(out_df.to_csv(index=False, header=False, sep="\n",
                                  quoting=csv.QUOTE_NONE, quotechar="'", escapechar="\\").replace("\\", ""))

    def to_percentile(self, nth, scope="global", inplace=False):
        self._to_dataframe()
        print(f"==> Transforming to {nth} percentile")
        ret_df = self.wiggle_df
        seqids = ret_df["variableStep_chrom"].unique().tolist()
        for seqid in seqids:
            col = ret_df[ret_df["variableStep_chrom"] == seqid]["score"]
            if self.orientation == "r":
                col = col.abs()
                if scope == "global":
                    ret_df.loc[ret_df["variableStep_chrom"] == seqid, "score"] = (col / np.percentile(col, nth)) * -1
                elif scope == "stretch":
                    pass
                    # TODO
                else:
                    print("Bad option")
            else:
                if scope == "global":
                    ret_df.loc[ret_df["variableStep_chrom"] == seqid, "score"] = (col / np.percentile(col, nth))
                elif scope == "stretch":
                    pass
                    # TODO
                else:
                    print("Bad option")

        ret_df["score"] = ret_df["score"].replace([np.nan, np.inf, -np.inf], 0.0)
        if inplace:
            self.wiggle_df = ret_df
            del ret_df
        else:
            return ret_df

    def to_step_height(self, step_range, step_direction, inplace=False):
        self._to_dataframe()
        print("==> Transforming to step height")
        ret_df = self.wiggle_df
        seqids = ret_df["variableStep_chrom"].unique().tolist()
        for seqid in seqids:
            if self.orientation == "r":
                ret_df.loc[ret_df["variableStep_chrom"] == seqid, "score"] = \
                    self._generate_step_height_col(ret_df[ret_df["variableStep_chrom"] == seqid]["score"].abs(),
                                                   step_range, step_direction, self.orientation) * -1
            else:
                ret_df.loc[ret_df["variableStep_chrom"] == seqid, "score"] = \
                    self._generate_step_height_col(ret_df[ret_df["variableStep_chrom"] == seqid]["score"],
                                                   step_range, step_direction, self.orientation)
        if inplace:
            self.wiggle_df = ret_df
            del ret_df
        else:
            return ret_df

    @staticmethod
    def _generate_step_height_col(in_col, step_range, step_direction, orientation):
        df = pd.DataFrame()
        df["scores"] = in_col
        df["mean_before"] = df["scores"]
        df["mean_after"] = df["scores"].shift(-(step_range + 1))
        df["mean_before"] = df["mean_before"].rolling(step_range).mean()
        df["mean_after"] = df["mean_after"].rolling(step_range).mean()
        if step_direction == "end_start" and orientation == "f":
            df["step_height"] = df["mean_before"] - df["mean_after"]
        elif step_direction == "end_start" and orientation == "r":
            df["step_height"] = df["mean_after"] - df["mean_before"]
        elif step_direction == "start_end" and orientation == "f":
            df["step_height"] = df["mean_after"] - df["mean_before"]
        elif step_direction == "start_end" and orientation == "r":
            df["step_height"] = df["mean_before"] - df["mean_after"]
        else:
            print("Error")
            exit(1)
        df[df["step_height"] < 0] = 0.0
        df["step_height"] = df["step_height"].shift(1)
        df.fillna(0.0, inplace=True)
        df.rename(columns={"step_height": in_col.name}, inplace=True)
        return df[in_col.name]

    def to_log2(self, inplace=False):
        self._to_dataframe()
        print("==> Transforming to Log2")
        ret_df = self.wiggle_df
        if self.orientation == "r":
            ret_df["score"] = np.log2(ret_df["score"].abs().replace([0, 0.0], np.nan)) * -1
        else:
            ret_df["score"] = np.log2(ret_df["score"].replace([0, 0.0], np.nan))
        ret_df["score"] = ret_df["score"].replace([np.nan, np.inf, -np.inf], 0.0)
        if inplace:
            self.wiggle_df = ret_df
            del ret_df
        else:
            return ret_df

    def to_log10(self, inplace=False):
        self._to_dataframe()
        print("==> Transforming to Log10")
        ret_df = self.wiggle_df
        if self.orientation == "r":
            ret_df["score"] = np.log10(ret_df["score"].abs().replace([0, 0.0], np.nan)) * -1
        else:
            ret_df["score"] = np.log10(ret_df["score"].replace([0, 0.0], np.nan))
        ret_df["score"] = ret_df["score"].replace([np.nan, np.inf, -np.inf], 0.0)
        if inplace:
            self.wiggle_df = ret_df
            del ret_df
        else:
            return ret_df

    def arithmethic(self, opt, value, inplace=False):
        self._to_dataframe()
        ret_df = self.wiggle_df
        if opt == "add":
            print(f"==> Adding {value} to coverage")
            #TODO
        elif opt == "sub":
            print(f"==> Subtracting {value} from coverage")
            #TODO
        elif opt == "mul":
            print(f"==> Multiplying coverage by {value}")
            if self.orientation == "r":
                ret_df["score"] = ret_df["score"].abs().multiply(value).multiply(-1)
            else:
                ret_df["score"] = ret_df["score"].multiply(value)
            ret_df["score"] = ret_df["score"].replace([np.nan, np.inf, -np.inf], 0.0)
        elif opt == "div":
            print(f"==> Dividing coverage by {value} ")
            if self.orientation == "r":
                ret_df["score"] = (ret_df["score"].abs() / value) * -1
            else:
                ret_df["score"] = ret_df["score"] / value
            ret_df["score"] = ret_df["score"].replace([np.nan, np.inf, -np.inf], 0.0)
        else:
            print("Bad option")
            exit(1)
        if inplace:
            self.wiggle_df = ret_df
            del ret_df
        else:
            return ret_df

    def split_wiggle(self, by, output_dir=None):
        self._to_dataframe()
        if by == "seqid":
            print(f"==> Splitting {os.path.basename(self.file_path)} by sequence ID")
            seqid_list = self.wiggle_df["variableStep_chrom"].unique()
            for seqid in seqid_list:
                if output_dir is None:
                    out_file_name = f"{os.path.dirname(self.file_path)}/{seqid}_" \
                                    f"{os.path.basename(self.file_path)}"
                else:
                    out_file_name = f"{os.path.abspath(output_dir)}/{seqid}_" \
                                    f"{os.path.basename(self.file_path)}"
                self.write_wiggle(out_file_name,
                                  alt_wiggle_df=self.wiggle_df[self.wiggle_df["variableStep_chrom"] == seqid])
        elif by == "fasta":
            print(f"==> Splitting {os.path.basename(self.file_path)} by fasta files")
            fasta_list = list(set([item["fasta"] for item in self.chrom_sizes]))
            for fasta in fasta_list:
                # get list of seqids for each fasta
                seqids = [x["seqid"] for x in self.chrom_sizes if fasta == x["fasta"]]
                if output_dir is None:
                    out_file_name = f"{os.path.dirname(self.file_path)}/{os.path.splitext(fasta)[0]}_"\
                                    f"{os.path.basename(self.file_path)}"
                else:
                    out_file_name = f"{os.path.abspath(output_dir)}/{os.path.splitext(fasta)[0]}_"\
                                    f"{os.path.basename(self.file_path)}"
                self.write_wiggle(out_file_name,
                                  alt_wiggle_df=self.wiggle_df[self.wiggle_df["variableStep_chrom"].isin(seqids)])
        else:
            print("Error: bad argument")
            exit(1)


"""
            for item in self.raw_data:
                if output_dir is None:
                    out_file_name = f"{os.path.dirname(self.file_path)}/{item['variableStep_chrom']}_"\
                                    f"{os.path.basename(self.file_path)}"
                else:
                    out_file_name = f"{os.path.abspath(output_dir)}/{item['variableStep_chrom']}_"\
                                    f"{os.path.basename(self.file_path)}"
                with open(out_file_name, 'w') as f:
                    print(f"==> Writing file: {os.path.basename(out_file_name)}")
                    f.write(f'track type={item["track_type"]} name="{item["track_name"]}"\n' + \
                            f'variableStep chrom={item["variableStep_chrom"]} span="{item["variableStep_span"]}"\n' + \
                            f'{nl.join([f"{k} {v}" for k, v in item["data"].items()])}')
"""