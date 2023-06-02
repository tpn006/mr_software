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
    
    def calcAndSetThreshold(self, background_sequences, pval):
        """
        Given a p-value threshold determine what score that is and set the threshold variable to that score
        """
        self.pvalThresh = pval
        for seq in background_sequences:
            self.scoreRandSequence(str(seq))
        #print(self.name + ' length of rnd array: ' + str(len(self.rnd_scr_arr)))
        #print('length of score array: ' + str(len(self.score_arr)))
        num = int(len(self.rnd_scr_arr) * (1-pval))
        self.rnd_scr_arr.sort()
        self.threshold_score = self.rnd_scr_arr[num]
        # print(self.rnd_scr_arr)
        self.setPval()

        '''i = 1
        temp_dist = self.rnd_scr_arr.copy()
        while i < num:
            temp_dist.pop(temp_dist.index(max(temp_dist)))
            i+=1
        thresh = max(temp_dist)
        self.threshold = thresh'''
        #self.graphThreshold()
        #print(str("SCORE ARRAY") + str(self.rnd_scr_arr))

    
    """    def setPval(self):
        peak_total = 0
        peak_motif = 0
        bg_total = 0

        bg_motif  = 0 #FIX THIS LATER
        # TODO figure out why the self.rnd_scr_arr is so long
        peak_total = len(self.score_arr)
        for score in self.score_arr:
            if score >= self.threshold:
                peak_motif += 1
        for score in self.rnd_scr_arr:
            if score >= self.threshold:
                bg_motif += 1
        #print("Running scipy with threshold "+ str(self.threshold) +" Peak total: " + str(peak_total) + " peak motif: " + str(peak_motif) + " Background total: " + str(peak_total) + " bgrnd motif: " + str(bg_motif))
        self.pval = KnownMotif.ComputeEnrichment(peak_total, peak_motif, peak_total, bg_motif)
    """
    def quickScore(self, sequence):
        if sequence == None: # IF we get grabage input :(
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

    def magicDoAllFunction(self, peak_seqs, back_seqs, p_val = 0.01):
        # Get the threshold score by scoring all background sequences
        back_scores = []
        for seq in back_seqs:
            back_scores.append(self.quickScore(seq))
        back_scores.sort()
        num_expected_to_pass_at_pval = int(len(back_scores)*(1-p_val))
        # Set the threshold score
        self.threshold_score = back_scores[num_expected_to_pass_at_pval]

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
        # table = [[num peaks that pass threshold, num peaks that don't pass],[num background that pass, num background that doesnt pass]]
        table = [[sig_seq, len(peak_scores)-sig_seq], [num_expected_to_pass_at_pval, len(back_scores)-num_expected_to_pass_at_pval]]
        print("Our fisher exact table for" + self.name+ " shows: " + str(table))
        odds, pval = scipy.stats.fisher_exact(table)
        self.pval = pval
        


        


    def setPval(self):
        #Assume threshold already set
        fwd_peak_pass = 0
        for seq in self.score_arr:
            if seq > self.threshold_score:
                fwd_peak_pass += 1  
        peak_total = len(self.score_arr)
        bg_total = peak_total # FOR NOW BECAUSE rnd_scr_arr is growing for no reason
        bg_motif= len(self.rnd_scr_arr) * self.pvalThresh
        table = [[fwd_peak_pass, peak_total-fwd_peak_pass], [bg_motif, bg_total-bg_motif]]
        odds, pval = scipy.stats.fisher_exact(table)
        self.pval = pval
        #self.pval = self.ComputeEnrichment(len(self.score_arr), fwd_peak_pass, len(self.score_arr), num_bg_pass)

    """
    def updateScore(self, val, which_score):
        
        #updates score by appending the value
        
        #print('before update score: ' + str(len(self.rnd_scr_arr)))
        if "forward" in which_score:
            self.score_arr.append(float(val))
        elif "random" in which_score:
            self.rnd_scr_arr.append(float(val))
        #print('after update score ' + str(len(self.rnd_scr_arr)))
    """

    def updateScore(self, val):
        self.score_arr.append(float(val))

    def updateScoreRand(self, val):
        self.rnd_scr_arr.append(float(val))

    def scoreSequence(self, sequence):
        """
        Given a long sequence from a peak, find where this motif is best suited and then update the score.
        Calls a lot of scoreSeq with shorter sequences
        """
        # list of scores. scores[i] should give the score of the substring sequence[i:i+n]. The default value for the score will be the maximum negative value so that the highest will always overwrite
        #scores = [0]*(len(sequence)-self.wscoreRandSequence + 1)
        scores = []
        i = 0
        # score the PWM on all spots in the sequence
        while i < ( len(sequence) - self.w + 1 ):
            scores.append(self.maxScore(sequence[i : i + self.w]))
            i += 1
        scores.sort()
        if len(scores) != 0:
            self.updateScore(scores[len(scores)-1]) #update score with the highest score found (last thing when sorted is the highest score)
        
    def scoreRandSequence(self, sequence):
        """
        Given a long sequence from a peak, find where this motif is best suited and then update the score.
        Calls a lot of scoreSeq with shorter sequences
        """
        # list of scores. scores[i] should give the score of the substring sequence[i:i+n]. The default value for the score will be the maximum negative value so that the highest will always overwrite
        #scores = [0]*(len(sequence)-self.w + 1)
        scores = []
        i = 0
        # score the PWM on all spots in the sequence
        while i < ( len(sequence) - self.w + 1 ):
            scores.append(self.maxScore(sequence[i : i + self.w -1]))
            i += 1
        scores.sort()
        if len(scores) != 0:
            self.updateScoreRand(max(scores)) #update score with the highest score found (last thing when sorted is the highest score)
    # WHY SITN THSI WORKING

    def ReverseComplement(self, sequence):
        """ Get the reverse complement of a sequence
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
        """
        while i >= 0:
            if sequence[i] == 'A':
                revcomp = revcomp + 'T'
            elif sequence[i] == 'C':
                revcomp = revcomp + 'G'
            elif sequence[i] == 'G':
                revcomp = revcomp + 'C'
            elif sequence[i] == 'T':
                revcomp = revcomp + 'A'
            elif sequence[i] == 'N':
                revcomp = revcomp + 'N'
            i = i-1
        return revcomp
        """

    def maxScore(self, sequence):
        """
        Score the sequence and its reverse complement and determine which is better
        """
        fwd_score = self.scoreSeq(sequence)
        bwd_score = self.scoreSeq(self.ReverseComplement(sequence))
        #print("Forward: " + str(fwd_score) + " Backward: " + str(bwd_score))
        if fwd_score >= bwd_score:
            return fwd_score
        else:
            return bwd_score

    def scoreSeq(self, sequence):
        """
        Given a sequence, find the score of that given sequence and return that score
        """
        score_temp = 0
        i = 0
        for nuc in str(sequence):
            # if the thing is an N, then the score is a 0
            if 'N' in nuc:
                score_temp += 0 
            else:
                print(i)
                print(nuc)
                score_temp += float(self.pwm[self.nucs[nuc]][i])
            i += 1
        return score_temp
    
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
        #self.setPval() # MOVE THIS LINE
        return "MOTIF: " + self.name + " threshold score: " + str(self.threshold_score) + " pval:" + str(self.pval)