#! /usr/bin/python3

import os


class BinaryBlast():
    '''
    A class for reading compiled BLAST databases.
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

    def __read_headers(self):
        """
        Reads header file, gets gi (s) from 'visible string' field and makes a dictionary with seq offset tuples
        """
        import re
        pattern=re.compile(b'gi\|(\d{1,9})\|')
        self.__header_dict = {}
        for j in range(len(self.__header_offsets)):
            self.headers.seek(self.__header_offsets[j][0])
            b=self.headers.read(self.__header_offsets[j][1])
            found=pattern.search(b)
            for a in found.groups():
                #There may be more than one gi for a given sequence. Nearly always so in nr, for example
                self.__header_dict[a] = self.__sequence_offsets[j]

    def __get_position(self, seqid):
        """
        Returns position of a seq with given seqid, either from dict or (way slower) scroll through headers
        :param seqid: int
        :return:
        """
        if self.load_headers:
            return self.__header_dict[seqid]
        else:
            for j in range(len(self.__header_offsets)):
                self.headers.seek(self.__header_offsets[j][0])
                b = self.headers.read(self.__header_offsets[j][1])
                if bytes('gi|'+str(seqid)+'|', 'ascii') in b:#Just bytes(str(seqid)) is a collision risk
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
        Returns a protein sequence of a given seqid as IUPAC string
        :param seqid:int
        :return:
        '''
        return self.__get_seq_by_position(self.__get_position(seqid))

    #Here will be iterator support (when I get to it)
    def __next__(self):
        pass


