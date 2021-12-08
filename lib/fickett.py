import numpy

class Fickett:
    '''
    calculate Fickett TESTCODE for full sequence
    NAR 10(17) 5303-531
    modified from source code of CPAT 1.2.1 downloaded from https://sourceforge.net/projects/rna-cpat/files/?source=navbar 
    '''
    def __init__(self):
        '''new compiled Fickett look-up table'''
        self.position_parameter  = [1.9,1.8,1.7,1.6,1.5,1.4,1.3,1.2,1.1,0.0]
        self.content_parameter  = [0.33,0.31,0.29,0.27,0.25,0.23,0.21,0.19,0.17,0]
        self.position_probability = {
            "A":[0.51,0.55,0.57,0.52,0.48,0.58,0.57,0.54,0.50,0.36],
            "C":[0.29,0.44,0.55,0.49,0.52,0.60,0.60,0.56,0.51,0.38],
            "G":[0.62,0.67,0.74,0.65,0.61,0.62,0.52,0.41,0.31,0.17],
            "T":[0.51,0.60,0.69,0.64,0.62,0.67,0.58,0.48,0.39,0.24],
            }
        self.position_weight = {"A":0.062,"C":0.093,"G":0.205,"T":0.154}
        self.content_probability = {
            "A":[0.40,0.55,0.58,0.58,0.52,0.48,0.45,0.45,0.38,0.19],
            "C":[0.50,0.63,0.59,0.50,0.46,0.45,0.47,0.56,0.59,0.33],
            "G":[0.21,0.40,0.47,0.50,0.52,0.56,0.57,0.52,0.44,0.23],
            "T":[0.30,0.49,0.56,0.53,0.48,0.48,0.52,0.57,0.60,0.51]
            }
        self.content_weight = {"A":0.084,"C":0.076,"G":0.081,"T":0.055}


    def look_up_position_probability(self,value, base):
        '''
        look up positional probability by base and value
        '''
        if float(value) < 0:
            return None
        for idx,val in enumerate (self.position_parameter):
            if (float(value) >= val):
                return float(self.position_probability[base][idx]) * float(self.position_weight[base])

    def look_up_content_probability(self,value, base):
        '''
        look up content probability by base and value
        '''
        if float(value) < 0:
            return None
        for idx,val in enumerate (self.content_parameter):
            if (float(value) >= val):
                return float(self.content_probability[base][idx]) * float(self.content_weight[base])

    def fickett_value(self,dna):
        '''
        calculate Fickett value from full RNA transcript sequence
        '''
        if len(dna) < 2:
            return 0
        fickett_score=0
        dna=dna
        total_base = len(dna)
        A_content = float(dna.count("A"))/total_base
        C_content = float(dna.count("C"))/total_base
        G_content = float(dna.count("G"))/total_base
        T_content = float(dna.count("T"))/total_base

        phase_0 = dna[::3]
        phase_1 = dna[1::3]
        phase_2 = dna[2::3]
        
        phase_0_A = phase_0.count("A")
        phase_1_A = phase_1.count("A")
        phase_2_A = phase_2.count("A")
        phase_0_C = phase_0.count("C")
        phase_1_C = phase_1.count("C")
        phase_2_C = phase_2.count("C")
        phase_0_G = phase_0.count("G")
        phase_1_G = phase_1.count("G")
        phase_2_G = phase_2.count("G")
        phase_0_T = phase_0.count("T")
        phase_1_T = phase_1.count("T")
        phase_2_T = phase_2.count("T")

        A_content = float(phase_0_A + phase_1_A + phase_2_A)/total_base
        C_content = float(phase_0_C + phase_1_C + phase_2_C)/total_base
        G_content = float(phase_0_G + phase_1_G + phase_2_G)/total_base
        T_content = float(phase_0_T + phase_1_T + phase_2_T)/total_base
        A_position= numpy.max([phase_0_A,phase_1_A,phase_2_A])/(numpy.min([phase_0_A,phase_1_A,phase_2_A]) +1.0)
        C_position= numpy.max([phase_0_C,phase_1_C,phase_2_C])/(numpy.min([phase_0_C,phase_1_C,phase_2_C]) +1.0)
        G_position= numpy.max([phase_0_G,phase_1_G,phase_2_G])/(numpy.min([phase_0_G,phase_1_G,phase_2_G]) +1.0)
        T_position= numpy.max([phase_0_T,phase_1_T,phase_2_T])/(numpy.min([phase_0_T,phase_1_T,phase_2_T]) +1.0)

        fickett_score += self.look_up_content_probability(A_content,"A")
        fickett_score += self.look_up_content_probability(C_content,"C")
        fickett_score += self.look_up_content_probability(G_content,"G")
        fickett_score += self.look_up_content_probability(T_content,"T")
        
        fickett_score += self.look_up_position_probability(A_position,"A")
        fickett_score += self.look_up_position_probability(C_position,"C")
        fickett_score += self.look_up_position_probability(G_position,"G")
        fickett_score += self.look_up_position_probability(T_position,"T")
            
        return fickett_score

