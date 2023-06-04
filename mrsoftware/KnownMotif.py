import sys
import scipy

class KnownMotif:
    score_arr = []
    pvalThresh = 0
    rnd_scr_arr = []
    threshold_score = 0 # The score that divides matching from not matching in null dist
    name = "" # Name of the known motif
    w = 0 # Length of the known motif
    alength = 4 # number of nucleotides in the pwm
    pwm = None # pwm[A C G T][pos] returns value at position of nucleotide 
    nucs = {"A": 0, "C": 1, "G": 2, "T": 3}
    pval = 0 # pval from fisher exact test, significance of enrichment

    def __init__(self, name, alength, w, pwm):
        self.name = name
        self.w = w
        self.alength = alength
        self.pwm = pwm
        self.score_arr = []
        self.rnd_scr_arr = []

    def quickScore(self, sequence):
        if sequence == None: # IF we get bad input :(
            return
        scores = []
        # Make a shifting window to go through the pwm with
        num_sequences = len(sequence) - self.w
        for shift in range(num_sequences):
            score = 0
            # Go through each letter in the sequence and add the score from that value in the PWM
            for i in range(self.w -1):
                nuc = str.upper((sequence[shift:shift + self.w])[i])
                if 'N' in nuc:
                    score += 0
                else:
                    score += self.pwm[self.nucs[nuc]][i]
            # Add the score to an array for each sequence
            scores.append(float(score))
        # Only return the max
        if scores == None or len(scores) < 1: #Lowkey idk why this is happening
            return
        return max(scores)

    
    def scoreAllGetThresholdScoreCalcPValue(self, peak_seqs, back_seqs, p_val = 0.01): 
        """
        Scores all sequences (peak seqs and random background seqs, determines threshold score, and calculates p value for motif)
        
        Parameters:
        -----------
        peak_seqs: array of strs
            list of all peak sequences from genome
        back_seqs: array of strs
            list of all randomly generated sequences
        p_val: float
            number used to calculate threshold score

        Returns:
        --------
        Nothing
        """
        
        # Get the threshold score by scoring all background sequences
        print("Evaluating " + str(self.name))
        back_scores = []
        for seq in back_seqs:
            back_scores.append(self.quickScore(seq))
        back_scores.sort()
        num_expected_to_fail_at_pval = int(len(back_scores)*(1-p_val))
        num_expected_to_pass_at_pval = int(len(back_scores)*(p_val))
        
        # Set the threshold score
        self.threshold_score = back_scores[num_expected_to_fail_at_pval]

        # score peak seqs
        sig_seq = 0 # This is the number of significant seqs
        peak_scores = []
        i = 0
        for seq in peak_seqs:
            # find the max of the forward and reverse
            if seq == None:
                continue
            fw_score = self.quickScore(seq)
            bw_score = self.quickScore(self.ReverseComplement(seq))
            if fw_score == None:
                print(seq + "has none type fw")
                continue
            if bw_score == None:
                print(seq + "has none type bw")
                continue
            curr_score = max([fw_score, bw_score ])
            peak_scores.append(curr_score)
            if curr_score >= self.threshold_score:
                sig_seq += 1

        # p value the stuff
        # table = [[num peaks that pass threshold, num peaks that don't pass],
        #          [num background that pass, num background that doesnt pass]]
        table = [[sig_seq, len(peak_scores)-sig_seq],
                 [num_expected_to_pass_at_pval, len(back_scores)-num_expected_to_pass_at_pval]]
        # print("Our fisher exact table for " + self.name+ " shows: " + str(table))
        odds, pval = scipy.stats.fisher_exact(table)
        self.pval = pval 
        print("Done evaluating " + self.name + " resulting pval is " + str(self.pval))

    def ReverseComplement(self, sequence):
        """
        Get the reverse complement of a sequence
                
        Parameters
        ----------
        sequence : str
            Sequence of nucleotides
            
        Returns
        -------
        revcomp : str
            Reverse complement of sequence
            """
        revcomp = ""
        sequence = str.upper(sequence)
        i = len(sequence)-1
        nuc_comps = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
        for nuc in reversed(sequence):
            if nuc in nuc_comps:
                revcomp += nuc_comps[str.upper(nuc)]
            else:
                revcomp += 'N'
                
        #print(str(len(sequence)) + " compared to " + str(len(revcomp)))
        return revcomp
    
    def printName(self):
        print(self.name)

    def printPWM(self):
        nucs = ["A","C","G","T"]
        for row in range(len(self.pwm)):
            print(nucs[row], end = ": ")
            for letter in self.pwm[row]:
                print(str(letter)[0:4], end = " ")
            print()

    def printScoreArr(self):
        print(str(self.score_arr))
    
    def getThresh(self):
        return self.threshold_score

    def getScoreArr(self):
        return self.score_arr

    def getPval(self):
        return self.pval
        
    # To String method    
    def __str__(self):
        return "MOTIF: " + self.name + " threshold score: " + str(self.threshold_score) + " pval:" + str(self.pval)