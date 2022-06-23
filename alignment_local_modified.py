import numpy as np

class AlignmentLocalModified:
    def __init__(self, seq1:str, seq2 : str):
        self.sub_matrix: dict[tuple[str, str], int]
        self._seq1: str = seq1
        self._seq2: str = seq2
        self.score_matrix = np.empty((seq2.__len__()+1, seq1.__len__()+1), dtype=int)
        self.tb_matrix = np.empty((seq2.__len__()+1, seq1.__len__()+1), dtype=tuple)
        self.gap_penalty: int = 1
        self.threshold_value: int = 10
        self.last: int
        self.last_tuple: tuple
        self.tb_seq1: str = seq1
        self.tb_seq2: list = list("-" * seq1.__len__())

    def set_gap(self, gap_penalty: int):
        self.gap_penalty = gap_penalty

    def set_threshold(self, threshold_value: int):
        self.threshold_value = threshold_value

    def align(self):
        if self.sub_matrix is None:
            raise self.MatrixUndefined("Substitution matrix hasn't been defined")
        else:
            self.initialize_score_matrix()
            previous_max = 0
            previous_max_tuple = (0, 0)
            for jcol in range(1, self._seq1.__len__()+1):
                current_max = 0
                current_max_tuple = (0, 0)
                for irow in range(0, self._seq2.__len__()+1):
                    if irow == 0:
                        self.score_matrix[irow, jcol] = max(previous_max - self.threshold_value,
                                                            self.score_matrix[irow, jcol - 1])

                        if self.score_matrix[irow, jcol] == previous_max - self.threshold_value:
                            self.tb_matrix[irow, jcol] = previous_max_tuple
                        else:
                            self.tb_matrix[irow, jcol] = (irow, jcol - 1)

                        if self.score_matrix[irow, jcol] > current_max:
                            current_max = self.score_matrix[irow, jcol]
                            current_max_tuple = (irow, jcol)
                    else:
                        this_row_max = self.score_matrix[0, jcol]
                        try:
                            diagonal_sc = self.score_matrix[irow - 1, jcol - 1] + self.sub_matrix[(self._seq1[jcol - 1],
                                                                                                   self._seq2[irow - 1])]
                        except KeyError:
                            diagonal_sc = self.score_matrix[irow - 1, jcol - 1] + self.sub_matrix[(self._seq2[irow - 1],
                                                                                                   self._seq1[jcol - 1])]
                        horizontal_sc = self.score_matrix[irow, jcol - 1] - self.gap_penalty
                        vertical_sc = self.score_matrix[irow - 1, jcol] - self.gap_penalty

#                        print(this_row_max, diagonal_sc, horizontal_sc, vertical_sc, sep='|')
                        self.score_matrix[irow, jcol] = max(this_row_max, diagonal_sc, horizontal_sc, vertical_sc)

                        if self.score_matrix[irow, jcol] == this_row_max:
                            self.tb_matrix[irow, jcol] = (0, jcol)
                        elif self.score_matrix[irow, jcol] == diagonal_sc:
                            self.tb_matrix[irow, jcol] = (irow - 1, jcol - 1)
                        elif self.score_matrix[irow, jcol] == horizontal_sc:
                            self.tb_matrix[irow, jcol] = (irow, jcol - 1)
                        elif self.score_matrix[irow, jcol] == vertical_sc:
                            self.tb_matrix[irow, jcol] = (irow - 1, jcol)

                        if self.score_matrix[irow, jcol] > current_max:
                            current_max = self.score_matrix[irow, jcol]
                            current_max_tuple = (irow, jcol)

                previous_max = current_max
                previous_max_tuple = current_max_tuple

            self.last = max(previous_max - self.threshold_value, self.score_matrix[0, self._seq1.__len__()])
            if self.last == previous_max - self.threshold_value:
                self.last_tuple = previous_max_tuple
            else:
                self.last_tuple = (0, self._seq1.__len__())

    def print_traceback_matrix(self):
        max_len = 0
        for row in self.tb_matrix:
            for el in row:
                el_in_str_len = str(el).__len__()
                if el_in_str_len > max_len:
                    max_len = el_in_str_len

        iter = -1
        print(end=' ' * (max_len + 4))
        for sym in self._seq1:
            print(sym, end='')
            print(' ' * (max_len - 1), end='|')
        print()

        for row in self.tb_matrix:
            if iter != -1:
                print(self._seq2[iter], end=' |')
            else:
                print(' ' * 3, end='')
            iter += 1
            for el in row:
                print(el, end='')
                print(' ' * (max_len - str(el).__len__()), end='|')
            print()

    def traceback(self):
        current_i, current_j = self.last_tuple
#        print(f'i: {current_i}, j:{current_j}')
        for iter in range(0, self._seq1.__len__()):
#            print(f'i: {current_i}, j:{current_j}, tb[i, j]: {self.tb_matrix[current_i, current_j][0]}')
            if (self.tb_matrix[current_i, current_j][0] != current_i) and (self.tb_matrix[current_i, current_j][0] != current_i - 1):
                self.tb_seq2[current_j-1] = '*'
#                print('*', end='')
            elif self.tb_matrix[current_i, current_j][0] == current_i or self.tb_matrix[current_i, current_j][1] == current_j:
                self.tb_seq2[current_j - 1] = '-'
#                print('-', end='')
            else:
                self.tb_seq2[current_j - 1] = self.tb_seq1[current_j - 1]
#                print(self.tb_seq1[current_j - 1], end='')

            current_i, current_j = self.tb_matrix[current_i, current_j]

    def print_traceback(self):
        print("Alignment of sequence1 and sequence2")
        print('\t',self.tb_seq1)
        print('\t',''.join(self.tb_seq2))

        irow = 0
        for jcol in range(self._seq1.__len__(), 0, -1):
            if jcol == self._seq1.__len__():
                pass


    def print_matrix(self):
        max_len = 0
        for row in self.score_matrix:
            for el in row:
                el_in_str_len = str(el).__len__()
                if el_in_str_len > max_len:
                    max_len = el_in_str_len

        iter = -1
        print(end=' ' * (max_len + 4))
        for sym in self._seq1:
            print(sym, end='')
            print(' ' * (max_len - 1), end='|')
        print()

        for row in self.score_matrix:
            if iter != -1:
                print(self._seq2[iter], end=' |')
            else:
                print(' ' * 3, end='')
            iter += 1
            for el in row:
                print(el, end='')
                print(' ' * (max_len - str(el).__len__()), end='|')
            print()

    def define_matrix(self, matrix: dict[tuple[str, str], int]):
        self.sub_matrix = matrix

    def initialize_score_matrix(self):
        self.score_matrix[:, 0] = 0


    class MatrixUndefined(Exception):
        pass