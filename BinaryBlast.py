#! /usr/bin/python3

import os

class BinaryBlast():
    '''
    A class for reading compiled protein BLAST databases.
    Take filename without extension on creation
    Return given sequence by BinaryBlast.getSeq(seqID) or act as iterator
    '''

    PROTEIN_CODE = ['-', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
                    'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y', 'Z', 'U', '*', 'O', 'J']

    def __init__(self, name, load_headers=False):
        '''
        :param name: str
        :param load_headers: bool
        :return: self: BinaryBlast
        Upon initialisation, this class reads index file and creates lists of header and sequence offsets.
        If load_headers is set to true, it also loads complete offset by ID dictionary. It improves speed, but is
        RAM-costly. Work only with protein DBs (.pin/phr/psq), support for DNA will be added later
        '''
        self.load_headers = load_headers
        self.filepath = os.path.abspath(name)
        indexFile = self.filepath + '.pin'
        sequenceFile = self.filepath + '.psq'
        headersFile = self.filepath + '.phr'
        self.index = open(indexFile, mode='rb')
        self.sequence = open(sequenceFile, mode='rb')
        self.headers = open(headersFile, mode='rb')
        # Add some exception handling here ASAP
        #Current code just passes through an os exception, but some file checking should be done
        self.__header_offsets = []
        self.__sequence_offsets = []
        self.index.seek(8)
        t = int.from_bytes(self.index.read(4), 'big')
        self.title = self.index.read(t).decode(encoding='ascii')  #Take title length and read title
        s = int.from_bytes(self.index.read(4),
                           'big')  #Timestamp length. Fuck the timestamp itself, no one ever needs it
        self.index.seek(t + s + 16)
        self.sequence_number = int.from_bytes(self.index.read(4), 'big')
        self.index.seek(t + s + 32)  #We skip residue count and max sequence length
        #Adding header offsets
        i = int.from_bytes(self.index.read(4), 'big')
        for j in range(self.sequence_number):  #There are N+1 ints, we read them N times in cycle and once before
            a = int.from_bytes(self.index.read(4), 'big')
            self.__header_offsets.append((i, a - i))
            i = a
        #Exactly the same for sequence offsets
        i = int.from_bytes(self.index.read(4), 'big')
        for j in range(self.sequence_number):
            a = int.from_bytes(self.index.read(4), 'big')
            self.__sequence_offsets.append((i, a - i-1))#Last seq byte is unncecesary: that's NUL byte separator
            i = a
        if self.load_headers:
            self.__read_headers()

    #Various internal functions unavailable to the user

    def __read_headers(self):
        """
        Reads header file, gets 'visible string' field(s) and makes a dictionary with seq offset tuples
        Also skips 'BL_ORD_ID' fields, but not others
        """
        self.__header_dict = {}
        for j in range(len(self.__header_offsets)):
            self.headers.seek(self.__header_offsets[j][0])
            b=self.headers.read(self.__header_offsets[j][1])
            #string_arr=b.split(sep=b'\x1A') #Every element *starts* with a VisibleString
            s=_vs_to_str(b)
            for a in s:
                if 'BL_ORD_ID' in a:
                    continue
                self.__header_dict[a]=self.__sequence_offsets[j]

    def __get_position(self, seqid):
        """
        Returns position of a seq with given seqid, either from dict or (slower) scroll through headers
        :param seqid: str
        :return:
        """
        if self.load_headers:
            return self.__header_dict[seqid]
        else:
            for j in range(len(self.__header_offsets)):
                self.headers.seek(self.__header_offsets[j][0])
                b = self.headers.read(self.__header_offsets[j][1])
                if bytes(seqid, encoding='ascii') in b: #We don't care
                    return self.__sequence_offsets[j]



    def __get_seq_by_position(self, pos):
        """
        Gets sequence by its position within .psq file. The only parameter is a tuple (offset,length) in bytes
        :param pos: tuple
        :return:
        """
        self.sequence.seek(pos[0])
        byte_seq = self.sequence.read(pos[1])
        return ''.join(self.PROTEIN_CODE[a] for a in byte_seq)

    def get_seq(self, seqid):
        '''
        Returns a protein sequence of a given seqid as SeqRecord object
        :param seqid:int
        :return:
        '''
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.Alphabet import IUPAC
        simple_seq=Seq(self.__get_seq_by_position(self.__get_position(seqid)), IUPAC.protein)
        seq_obj = SeqRecord(simple_seq)
        seq_obj.id=seqid
        return seq_obj

    #Here will be iterator support (when I get to it)
    def __iter__(self):
        self.__iter_coord=0
        return self

    def __next__(self):
        if self.__iter_coord >= self.sequence_number:
            raise StopIteration
        self.headers.seek(self.__header_offsets[self.__iter_coord][0])
        h=self.headers.read(self.__header_offsets[self.__iter_coord][1])


def _vs_to_str(s):
    '''
    Returns all VisibleString elements in a given header as strings
    :param s: bytes
    :return:
    '''
    looked=0
    r=[]
    while True:
        start=s.find(b'\x1A', looked)
        if start==-1:
            break
        if not s[start+1]>>7:
            len_len=1
            str_len=s[start+1]
            vs=s[start+2:start+len_len+str_len+1].decode(encoding='ascii')
        else:
            len_len=s[start+1]&0b01111111
            str_len=int.from_bytes(s[start+2:start+len_len+2], 'big')
            try:
                vs=s[start+len_len+2:start+len_len+str_len+2].decode(encoding='ascii')
            except UnicodeDecodeError:
                print('Line:\n{0}\nis broken. Length header {1} bytes, presumed\
                      string length {2}, entire line length {3}. Start at \
                      {4}'.format(s, len_len, str_len, len(s), start))
                quit()
        looked=start+len_len+str_len+2
        r.append(vs)
    return r



