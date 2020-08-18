import pandas as pd
import os
import multiprocessing as mp


class WiggleMatrix:
    def __init__(self, parsed_wiggles, chrom_sizes, processes=None):
        self.parsed_wiggles = parsed_wiggles
        self.processes = processes
        self.chrom_sizes = chrom_sizes
        backbone = []
        for chrom in self.chrom_sizes:
            backbone = [[chrom["seqid"], loc] for loc in range(1, chrom["size"] + 1, 1)]
        self.wiggle_matrix_df = pd.DataFrame(backbone, columns=["seqid", "location"])
        self.build_matrix()
        self.f_wiggle_matrix_df = None
        self.r_wiggle_matrix_df = None

    def build_matrix(self):
        pool = mp.Pool(processes=self.processes)
        processes = []
        print("Building the matrix from the parsed files")
        for parsed_wiggle in self.parsed_wiggles:
            processes.append(pool.apply_async(self._merge_single_wiggle_to_matrix, (parsed_wiggle,)))
        columns_series = [p.get() for p in processes]
        for column in columns_series:
            self.wiggle_matrix_df[column.name] = column
        self.get_matrix_by_orientation()

    def _merge_single_wiggle_to_matrix(self, wig):
        condition_name = wig.at[0, "track_name"]

        self.wiggle_matrix_df = pd.merge(how='outer',
                                         left=self.wiggle_matrix_df,
                                         right=wig[["variableStep_chrom", "location", "score"]],
                                         left_on=['seqid', 'location'],
                                         right_on=['variableStep_chrom', 'location']).fillna(0.0)
        self.wiggle_matrix_df.rename(columns={"score": condition_name}, inplace=True)
        for column in self.wiggle_matrix_df.columns.tolist():
            if "variableStep_chrom" in column:
                self.wiggle_matrix_df.drop(column, axis=1, inplace=True)
        print(f"==> Merged condition {condition_name} to matrix")
        return self.wiggle_matrix_df[condition_name]

    def get_matrix_by_orientation(self):
        f_column_list = ["seqid", "location"]
        r_column_list = ["seqid", "location"]
        for column in self.wiggle_matrix_df.columns.tolist():
            if "seqid" != column != "location":
                if self.wiggle_matrix_df[self.wiggle_matrix_df[column] < 0].empty:
                    f_column_list.append(column)
                if self.wiggle_matrix_df[self.wiggle_matrix_df[column] > 0].empty:
                    r_column_list.append(column)
        self.f_wiggle_matrix_df = self.wiggle_matrix_df.loc[:, f_column_list]
        self.r_wiggle_matrix_df = self.wiggle_matrix_df.loc[:, r_column_list]

    def agg_merge(self, by=None):
        self.get_matrix_by_orientation()
        print(f"==> Aggregating by {by}")
        f_cond_cols = self.f_wiggle_matrix_df.columns.tolist()
        f_cond_cols = [x for x in f_cond_cols if x not in ["seqid", "location"]]
        r_cond_cols = self.r_wiggle_matrix_df.columns.tolist()
        r_cond_cols = [x for x in r_cond_cols if x not in ["seqid", "location"]]
        if by == "max":
            self.f_wiggle_matrix_df["agg_col_forward"] = self.f_wiggle_matrix_df.loc[:, f_cond_cols].max(axis=1)
            self.r_wiggle_matrix_df["agg_col_reverse"] = self.r_wiggle_matrix_df.loc[:, r_cond_cols].min(axis=1)
        elif by == "min":
            self.f_wiggle_matrix_df["agg_col_forward"] = self.f_wiggle_matrix_df.loc[:, f_cond_cols].min(axis=1)
            self.r_wiggle_matrix_df["agg_col_reverse"] = self.r_wiggle_matrix_df.loc[:, r_cond_cols].max(axis=1)
        elif by == "average":
            self.f_wiggle_matrix_df["agg_col_forward"] = self.f_wiggle_matrix_df.loc[:, f_cond_cols].mean(axis=1)
            self.r_wiggle_matrix_df["agg_col_reverse"] = self.r_wiggle_matrix_df.loc[:, r_cond_cols].mean(axis=1)
        else:
            print("Fatal error")

    @staticmethod
    def write_matrix_to_wiggle_files(matrix, out_dir, prefix=None):
        print("==> Writing wiggle files")
        seqids = matrix["seqid"].unique()
        columns = [col for col in matrix.columns if col not in ["seqid", "location"]]
        for col in columns:
            out_str = ""
            out_str += f'track type=wiggle_0 name="{col}"\n'
            for seqid in seqids:
                out_str += f'variableStep chrom={seqid} span=1\n'
                out_str += matrix[matrix["seqid"] == seqid][["location", col]]\
                    .to_csv(index=False, header=False, mode='a', sep=" ")
            out_dir = os.path.abspath(out_dir)
            file = open(os.path.abspath(f"{os.path.dirname(out_dir)}/{prefix}.wig"), "w")
            file.write(out_str)
            file.close()
            s = '\n     └── '
            print(f"===> Wrote file: {prefix}.wig, contains sequence IDs:\n"
                  f"     └── {s.join(seqids)}")
