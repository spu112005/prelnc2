class FindCDS:
    '''
    Find the most like CDS in a given sequence 
    The most like CDS is the longest ORF found in the sequence
    When having same length, the upstream ORF is printed
    modified from source code of CPAT 1.2.1 downloaded from https://sourceforge.net/projects/rna-cpat/files/?source=navbar
    '''
    def __init__(self,seq):
        self.seq = seq
        self.result = (0,0,0,0,0)
        self.longest = 0
        self.basepair = {"A":"T","T":"A","U":"A","C":"G","G":"C","N":"N","X":"X"}

    def _reversecompliment(self):
        return "".join(self.basepair[base] for base in self.seq)[::-1]

    def get_codons(self,frame_number):
        '''
        Record every nucleotide triplet and its coordinate position for input sequence in one frame
        '''
        coordinate = frame_number
        while coordinate + 3 <= len(self.seq):
            yield (self.seq[coordinate:coordinate+3], coordinate)
            coordinate += 3 
    
    def find_longest_in_one(self,myframe,direction,start_codon,stop_codon):
        '''
        find the longest ORF in one reading myframe
        '''
        triplet_got = self.get_codons(myframe)  
        starts = start_codon
        stops = stop_codon
        '''
        Extend sequence by triplet after start codon encountered
        End ORF extension when stop codon encountered
        '''
        while True:
            try: 
                codon,index = next(triplet_got)
            except StopIteration:
                break 
            if codon in starts and codon not in stops:
                '''
                find the ORF start
                '''
                orf_start = index
                end_extension = False
                while True:
                    try: 
                        codon,index = next(triplet_got)
                    except StopIteration:
                        end_extension = True
                        integrity = -1
                    if codon in stops:
                        integrity = 1
                        end_extension = True
                    if end_extension:
                        orf_end = index + 3
                        Length = (orf_end - orf_start)
                        if Length > self.longest:
                            self.longest = Length
                            self.result = [direction,orf_start,orf_end,Length,integrity]
                        if Length == self.longest and orf_start < self.result[1]:
                            '''
                            if ORFs have same length, return the one that if upstream
                            '''
                            self.result = [direction,orf_start,orf_end,Length,integrity]
                        break

    def longest_orf(self,direction,start_codon={"ATG":None}, stop_codon={"TAG":None,"TAA":None,"TGA":None}):
        return_orf = ""
        for frame in range(3):
            self.find_longest_in_one(frame,"+",start_codon,stop_codon)
        return_orf = self.seq[self.result[1]:self.result[2]][:]
        start_coordinate = self.result[1]
        strand_direction = "+"
        orf_integrity = self.result[4]
        '''
        Also check reverse chain if -r is chosen
        '''
        if direction == "-":
            self.seq = self._reversecompliment()
            for frame in range(3):
                self.find_longest_in_one(frame,"-",start_codon,stop_codon)
            if self.result[0] == "-":
                return_orf = self.seq[self.result[1]:self.result[2]][:]
                start_coordinate = self.result[1]
                strand_direction = "-"
                orf_integrity = self.result[4]
        return return_orf,start_coordinate,strand_direction,orf_integrity

