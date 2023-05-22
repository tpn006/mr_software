class KnownMotif:
    value = 0 # Running total of the highest score of each peak sequence
    name = "" # Name of the known motif
    w = 0 # Length of the known motif
    alength = 4 # number of nucleotides in the pwm
    pwm = [[0]*w]*alength # pwm[pos][A C G T] returns value at position of nucleotide
    nucs = {"A": 0, "C": 1, "G": 2, "T": 3}
    num_seqs_scored = 0

    def __init__(self, name, alength, w, pwm):
        self.name = name
        self.w = w
        self.alength = alength
        self.pwm = pwm
    
    def updateScore(self, val):
        self.value += val
        self.num_seqs_scored += 1

    def scoreSequence(self, sequence):
        """
        Given a long sequence from a peak, find where this motif is best suited and then update the score
        """
        scores = [0]*(len(sequence)-self.w + 1) # list of scores. scores[i] should give the score of the substring sequence[i:i+n]
        i = 0
        # score the PWM on all spots in the sequence
        while i < ( len(sequence) - self.w + 1 ):
            scores[i] = self.scoreSeq(self, sequence[i:i+self.w])
            i += 1
            
        self.updateScore(scores.sort()[:-1]) #update score with the highest score found (last thing when sorted is the highest score)
    
    def scoreSeq(self, sequence):
        """
        Given a sequence of the same length as the PWM, find the score of that given sequence
        """
        score = 0
        i = 0
        for nuc in sequence:
            score += self.pwm(i, self.nucs[nuc])
            i += 1
        return score