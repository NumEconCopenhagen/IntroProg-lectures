import numpy as np
from types import SimpleNamespace
from itertools import permutations




class Person:
    def __init__(self, id,u, ranking):
        '''
        Initiate person class
        '''

        self.id = id
        self.u = u # Utility of working with the partner, ordered by their id
        self.ranking = ranking  # Ordered preferred list of partner ids
        self.partner = None # Chosen partner
        self.k = 0 # Index of the preferred list the person has gotten to if they are the 
        self.k_max = len(ranking) # Maximum index of the preferred list

    def next_offer(self):
        '''
        Returns the id of the next person that the person wants to work with
        or None if the person has already found a partner or has gone through everybody
        '''

        if self.k == self.k_max or self.partner is not None:
            return None        
        else:
            self.k += 1
            return self.ranking[self.k-1]


    def react_to_offer(self,offers,printit = False):
        '''
        Offers is a list of offers from proposers to recievers 
        This functions goes through that list and finds all offers to self,
        then chooses the best offer (Others should be rejected)
        '''

        if self.partner is None:
            u_partner = 0
        else:
            u_partner = self.u[self.partner]

        for i,o in enumerate(offers):

            if o == self.id:
                if printit:
                    print(f'Offer from {i+1} is {o+1}')
                    print(f'Utility: {self.u[i]:.4f} vs {u_partner:.4f}')

                if self.u[i]>u_partner:
                    if printit:
                        print(f'Accepting offer from {i+1} so far')
                    # If the offer is better than the current partner, then choose the offer
                    self.partner = i
                    u_partner = self.u[i]      

        return self.partner

    def __repr__(self):
        '''
        Text representation
        '''
        if self.partner is None:
            return f'Person {self.id+1} works alone'
        else:
            return f'Person {self.id+1} works with {self.partner+1}'

class MatchingModel:

    def __init__(self,S=10,M=10):
        '''
        Initiate the matching model
        '''
        par = self.par = SimpleNamespace()
        sol = self.sol = SimpleNamespace()
        par.S = S
        par.M = M

    def simulate_preferences(self):
        '''
        Simulate preferences for students and mentors
        '''
        
        par = self.par
        S = par.S
        M = par.M

        # Simulate preferences
        par.S_pref = np.random.uniform(size=(S,M))
        par.M_pref = np.random.uniform(size=(M,S))
        
        # Create ranking
        # In S_ranking each row is a student, and the columns are their rankings ordered from left to right.
        # So having 3 in the first column means prefering the fourth mentor the most
        
        par.S_ranking = np.argsort(-par.S_pref,axis=1)
        par.M_ranking = np.argsort(-par.M_pref,axis=1)

        # Create persons
        par.S_list = [Person(i, par.S_pref[i,:], par.S_ranking[i,:]) for i in range(S)]
        par.M_list = [Person(i, par.M_pref[i,:], par.M_ranking[i,:]) for i in range(M)]

    def reset(self):
        '''
        Reset person counter for everybody
        '''        

        par = self.par
        for p in par.S_list+par.M_list:
            p.k = 0
            p.partner = None
    
    
    
    def print_matching(self):

        par = self.par
        print('Matching')

        print('Students:')
        for p in par.S_list:
            print(p)

        print('\nMentors')
        for p in par.M_list:
            print(p)

        matchprint= '{'
        for x,y in self.current_matching():
            if x is not None:
                x += 1
            if y is not None:   
                y += 1
            matchprint += f'({x},{y}), '
        matchprint += '}'
        
        print('{(s:m)} :',matchprint)

    def current_matching(self):
        '''
        Returns current matching in the form of a list of tuples (s,m) where s is the student and m is the mentor they are matched with
        Order by the index of the smallest group
        '''
        par = self.par
        
        if par.S < par.M:
            current_matching = [(s.id,s.partner) for s in par.S_list]
        else:
            current_matching = [(m.partner,m.id) for m in par.M_list]

        return current_matching

    def check_matching(self):
        '''
        Checks if the current matching is consistent, in the sense that everybody is matched with who they think that they are matched with,
         and whether the matching is stable by looking for blocking pairs
        '''

        par = self.par

        for s in par.S_list:
            
            if s.partner is None:
                continue
            else:
                m = par.M_list[s.partner]
                if m.partner != s.id:
                    print(f'Error: {s+1} is not matched with {m+1}')

        for m in par.M_list:

            if m.partner is None:
                continue
            else:
                s = par.S_list[m.partner]
                if s.partner != m.id:
                    print(f'Error: {m+1} is not matched with {s+1}')

        
        self.check_stability()


    def check_stability(self,print_blocking_pairs=True):
        '''
        
        Check for any blocking pairs in the current matching

        '''

        par = self.par
        for s in par.S_list:

            if s.partner is None:
                partner_u = 0
            else:
                partner_u = s.u[s.partner]
            
            for i in range(s.k_max):

                if s.u[i]>partner_u: # Does s prefer someone else more than their partner
                
                    m_other= par.M_list[i]
                
                    # Does this person also prefer s over their partner or don't have a partner
                    if not (m_other.partner is None):
                        m_other_u = m_other.u[m_other.partner]
                        if m_other.u[s.id] <= m_other.u[m_other.partner]: 
                            continue # No blocking pair

                    if print_blocking_pairs:
                        print(f'Error: student {s.id+1} and mentor {i+1} is a blocking pair')
                    return False
   
        return True

    def DAA(self, proposers='S',print_matching=True):
        '''
        Implement deferred acceptance algorithm
        '''

        par = self.par
        sol = self.sol 

        recievers = 'M' if proposers == 'S' else 'S'
        P_list = getattr(par,f'{proposers}_list' )
        R_list = getattr(par,f'{recievers}_list')

        # Reset
        self.reset()

        round = 1 
        while True:

            if round==1:
                print('Starting DA algorithm')
            else:
                print(f'Round {round}')
            


            offers = [p.next_offer() for p in P_list]
            
            if all([o is None for o in offers]):
                print('All proposers have an outstanding offer')
                break

            chosen_list = [r.react_to_offer(offers) for r in R_list]

            # Reset partners for proposers
            for p in P_list:
                p.partner = None

            # set the accepted offerers as having a partner (outstanding offer untill the algorithm is done )            
            for i, chosen in enumerate(chosen_list):
                if chosen is not None:
                    P_list[chosen].partner = i

            if print_matching:
                print('Current matching:', self.current_matching())
            
    
            round += 1
            if round >= 100:
                print('Breaking, more than a 100 rounds')
                break

        
        setattr(sol,proposers+'_DAA', self.current_matching())

        assert self.check_stability() ,'Something went wrong DAA matching is not stable'

        if print_matching:
            self.print_matching()
    

    def all_matchings(self):
        '''
        Find all possible matches between students and mentors
        '''

        par = self.par
        S = par.S
        M = par.M
        matchings = []


        # Find the biggest group (not important when S=M)
        if S >= M:
            x = list(range(S))
            y = list(range(M))
            reverse = False
        else:
            x = list(range(M))
            y = list(range(S))
            reverse = True

        # Find a list of all combinations
        for comb in permutations(x,len(y)):

            matching = list(zip(comb,y))
            if reverse:
                matching = [(s,m) for m,s in matching]

            matchings.append(matching)

        return matchings


    def find_all_stable_matches(self):
        '''
        Among all possible matches, find all that are stable
        '''

        par = self.par
        sol = self.sol

        S = par.S
        M = par.M

        all_matchings = self.all_matchings()
        sol.stable_matchings = []
          
        
        for matching in all_matchings:
            self.reset()
            for s,m in matching:
                par.S_list[s].partner = m
                par.M_list[m].partner = s
            if self.check_stability(print_blocking_pairs=False):
                
                matchprint= '{'
                for x,y in matching:
                    if x is not None:
                        x += 1
                    if y is not None:   
                        y += 1
                    matchprint += f'({x},{y}), '
                matchprint = matchprint[:-2]
                matchprint += '}'
                print(matchprint)

                sol.stable_matchings.append(matching)

    
    def calculate_utility(self):
        '''
        Calculate the average utility of each group in all stable matchings
        '''
        par = self.par
        sol = self.sol
        
        
        for matching in sol.stable_matchings:
            S_utility = 0
            M_utility = 0

            for s,m in matching:
                S_utility += par.S_pref[s,m]/par.S
                M_utility += par.M_pref[m,s]/par.M

            if matching == sol.S_DAA:
                print('DAA, S proposes:',f'-> S: {S_utility:.4f}, M: {M_utility:.4f}')
            elif matching == sol.M_DAA:
                print('DAA, M proposes:',f'-> S: {S_utility:.4f}, M: {M_utility:.4f}')
            else:
                matchprint= '{'
                for x,y in matching:
                    if x is not None:
                        x += 1
                    if y is not None:   
                        y += 1
                    matchprint += f'({x},{y}), '
                matchprint = matchprint[:-2]
                matchprint += '}'
                print(matchprint,f'-> S: {S_utility:.4f}, M: {M_utility:.4f}')
        