from wiggletools.wiggle import Wiggle
import pandas as pd
import os
from Bio import SeqIO
import multiprocessing as mp


class WiggleMatrix:
    def __init__(self, fasta_file, wiggle_files, processes=None):
        self.fasta_file = fasta_file
        self.wiggle_files = wiggle_files
        self.wiggle_matrix_df = None
        if processes is None:
            self.processes = mp.cpu_count()
        elif processes.isnumeric():
            processes = int(processes)
            if processes > mp.cpu_count():
                self.processes = mp.cpu_count()
            elif processes > 0:
                self.processes = processes
            else:
                print("Error in processes number")
                exit(1)
        else:
            print("Error in processes number")
            exit(1)


    def build_matrix(self):
        parsed_fasta = SeqIO.parse(self.fasta_file, "fasta")
        pool = mp.Pool(processes=self.processes)
        processes = []
        for wig in self.wiggle_files:
            processes.append(pool.apply_async(self._parse_single_wiggle, (wig, )))
        parsed_wiggles = [p.get() for p in processes]
        pool.close()
        base_columns_lst = []
        # Initialize matrix with full length genome
        total_len = 0
        # Build dataframe backbone
        for seq_record in parsed_fasta:
            base_columns_lst = [[seq_record.id, i] for i in range(1, len(seq_record.seq) + 1, 1)]
            total_len += len(seq_record.seq)
        self.wiggle_matrix_df = pd.DataFrame(data=base_columns_lst, columns=["seqid", "location"])
        pool = mp.Pool(processes=self.processes)
        processes = []
        print("Building the matrix from the parsed files")
        for parsed_wiggle in parsed_wiggles:
            processes.append(pool.apply_async(self._merge_single_wiggle_to_matrix, (parsed_wiggle,)))
        columns_series = [p.get() for p in processes]
        for column in columns_series:
            self.wiggle_matrix_df[column.name] = column
        return self.wiggle_matrix_df

    def get_matrix_by_orientation(self):
        self.build_matrix()
        f_column_list = ["seqid", "location"]
        r_column_list = ["seqid", "location"]
        for column in self.wiggle_matrix_df.columns.tolist():
            if "seqid" != column != "location":
                if self.wiggle_matrix_df[self.wiggle_matrix_df[column] < 0].empty:
                    f_column_list.append(column)
                if self.wiggle_matrix_df[self.wiggle_matrix_df[column] > 0].empty:
                    r_column_list.append(column)
        return self.wiggle_matrix_df[f_column_list], self.wiggle_matrix_df[r_column_list]

    @staticmethod
    def _parse_single_wiggle(wig):
        return Wiggle(wig).parse()

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
        return self.wiggle_matrix_df[condition_name]

    def write_matrix_to_wiggle_files(self, matrix, out_dir, prefix=None):
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
            file = open(os.path.abspath(f"{os.path.dirname(out_dir)}/{prefix}{col}.wig"), "w")
            file.write(out_str)
            file.close()
            s = '\n     └── '
            print(f"===> Wrote file: {prefix}{col}.wig, contains sequence IDs:\n"
              f"     └── {s.join(seqids)}")

    def to_percentile(self, nth):
        pass

    def to_step_height(self, step_range):
        f_wiggles_matrix, r_wiggles_matrix = self.get_matrix_by_orientation()
        f_pool = mp.Pool(processes=8)
        r_pool = mp.Pool(processes=8)
        f_processes = []
        r_processes = []
        df_cols = ["seqid", "location"]
        seqid_list = self.wiggle_matrix_df.seqid.unique()
        for seqid in seqid_list:
            for col in f_wiggles_matrix.columns:
                if col not in df_cols:
                    f_processes.append(f_pool.apply_async(
                        self._generate_step_factor_col,
                        (f_wiggles_matrix[f_wiggles_matrix.seqid == seqid][col], step_range, "f",)))
            for col in r_wiggles_matrix.columns:
                if col not in df_cols:
                    r_processes.append(r_pool.apply_async(
                        self._generate_step_factor_col,
                        (r_wiggles_matrix[r_wiggles_matrix.seqid == seqid][col].abs(), step_range, "r",)))
        f_done_process = [p.get() for p in f_processes]
        r_done_process = [p.get() for p in r_processes]
        f_cols = []
        r_cols = []
        for col in f_done_process:
            f_wiggles_matrix[col.name] = col
            f_cols.append(col.name)
        for col in r_done_process:
            r_wiggles_matrix[col.name] = col * -1
            r_cols.append(col.name)


    def to_log2(self):
        pass

    def to_log10(self):
        pass

    @staticmethod
    def _generate_step_factor_col(in_col, step_range, orientation):
        df = pd.DataFrame()
        df["scores"] = in_col
        df["mean_before"] = df["scores"]
        df["mean_after"] = df["scores"].shift(-(step_range + 1))
        df["mean_before"] = df["mean_before"].rolling(step_range).mean()
        df["mean_after"] = df["mean_after"].rolling(step_range).mean()
        if orientation == "f":
            df["step_factor"] = df["mean_before"] - df["mean_after"]
        elif orientation == "r":
            df["step_factor"] = df["mean_after"] - df["mean_before"]
        else:
            print("Error")
            exit(1)
        df[df["step_factor"] < 0] = 0.0
        df["step_factor"] = df["step_factor"].shift(1)
        df.fillna(0.0, inplace=True)
        df.rename(columns={"step_factor": in_col.name}, inplace=True)
        return df[in_col.name]