class Sequence:
    """
    a class for manipulating sequences
    """

    def __init__(self, seq_str, name=b"", quality_str=b""):
        """
        initing by help of sequence string and name and quality string
        """
        self.seq_str = seq_str
        self.size = len(seq_str)
        self.name = name
        self.quality_str = quality_str

    def get_name(self):
        """
        :return: name of sequence in ascii 
        """
        return self.name.decode("ascii")

    @staticmethod
    def read_sequence(file_name):
        """
        read a fasta file containing one or more sequences
        :return: a list of sequences
        """
        f = open(file_name)
        stream = f.read().encode("ascii")
        f.close()
        reads = stream.split(b">")
        l = [Sequence("")] * 0
        for read_i in range(1 if stream[0] == 62 else 0, len(reads)):
            read = reads[read_i]
            first_new_line = read.find(b"\n") if stream[0] == 62 else 0
            name = read[:first_new_line].strip()
            seq_str = read[first_new_line + 1:].replace(b"\n", b"").replace(b"A", b"\x00").replace(
                b"C", b"\x01").replace(b"G", b"\x02").replace(b"T", b"\x03").replace(b"N", b"\x04")
            l.append(Sequence(seq_str, name))
        return l

    def __str__(self):
        return "%s: %s" % (self.name, self.seq_str)

    def reverse(self):
        """
        :return: reverse of this sequence 
        """
        return Sequence(self.seq_str[::-1].replace(b"\x00", b"\x10").replace(b"\x03", b"\x00").replace(b"\x10", b"\x03")
                        .replace(b"\x01", b"\x20").replace(b"\x02", b"\x01").replace(b"\x20", b"\x02"))

    def sub_genome(self, start_pos, end_pos):
        """
        a genome that its sequence is a substring of this sequence
        :return: result sequence
        """
        return Sequence(self.seq_str[start_pos:end_pos], self.name, self.quality_str)

    def get_readable_str(self):
        """
        convert sequence string to readable format (from binary)
        :return: string
        """
        return self.seq_str.replace(b"\x00", b"A").replace(b"\x01", b"C").replace(b"\x02", b"G"). \
            replace(b"\x03", b"T").replace(b"\x04", b"N").replace(b"\x05", b"U").replace(b"\x06", b"R"). \
            replace(b"\x07", b"Y").replace(b"\x08", b"M").replace(b"\x09", b"K").replace(b"\x0A", b"W"). \
            replace(b"\x0B", b"S").replace(b"\x0C", b"B").replace(b"\x0D", b"D").replace(b"\x0E", b"H"). \
            replace(b"\x0F", b"V").decode("ascii")

    def local_global_align(self, sequence, pos=None, startswith=b""):
        """
        align all of a given sequence to local part of self. maybe defined by pos
        :return: a tuple of score, start position and end position
        """
        s1 = self.seq_str if not pos else self.seq_str[pos[0]:pos[1]]
        s2 = sequence.seq_str
        s = [[-i for i in range(len(s2) + 1)] for j in range(len(s1) + 1)]
        sp = [[j for i in range(len(s2) + 1)] for j in range(len(s1) + 1)]
        for i in range(len(s1)):
            for j in range(len(s2)):
                if s1[i] == s2[j]:
                    maxx = s[i][j] + 1
                    maxxp = sp[i][j]
                else:
                    maxx = s[i][j] - 1
                    maxxp = sp[i][j]
                if s[i + 1][j] - 1 > maxx:
                    maxx = s[i + 1][j] - 1
                    maxxp = sp[i + 1][j]
                if s[i][j + 1] - 1 > maxx:
                    maxx = s[i][j + 1] - 1
                    maxxp = sp[i][j + 1]
                s[i + 1][j + 1] = maxx
                sp[i + 1][j + 1] = maxxp
        maxx = -1000000000
        maxi = 0
        # check all of end places to find best place of endings
        len_startswith = len(startswith)
        for i in range(len(s1) + 1):
            spi = sp[i][len(s2)]
            if s[i][len(s2)] > maxx and s1[spi - len_startswith:spi] == startswith:
                maxx = s[i][len(s2)]
                maxi = i
        return maxx, pos[0] + sp[maxi][len(s2)], pos[0] + maxi

    def global_align(self, sequence, match_score=1):
        """
        global align a sequence by self. maybe receiving a match score
        :return: score of match
        """
        s1 = self.seq_str
        s2 = sequence.seq_str
        s = [[-max(i, j) for i in range(len(s2) + 1)] for j in range(len(s1) + 1)]
        for i in range(len(s1)):
            for j in range(len(s2)):
                if s1[i] == s2[j]:
                    maxx = s[i][j] + match_score
                else:
                    maxx = s[i][j] - 1
                if s[i + 1][j] - 1 > maxx:
                    maxx = s[i + 1][j] - 1
                if s[i][j + 1] - 1 > maxx:
                    maxx = s[i][j + 1] - 1
                s[i + 1][j + 1] = maxx
        return s[len(s1)][len(s2)]
