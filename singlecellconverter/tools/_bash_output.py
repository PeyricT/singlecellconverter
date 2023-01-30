class ProgressBar:

    def __init__(self, steps, width=50):
        self.steps = steps
        self.progress = 0.0
        self.max = max(steps.values())
        self.width = width
        self.iter = None

    def start(self):
        self.iter = iter(self.steps)
        bar = '|'+' '*self.width+'|'
        print(f'\r  {bar}    0% : starting...', end='\r')

    def update(self):
        try:
            key = next(self.iter)
        except StopIteration:
            return

        self.progress = float(self.steps[key]/self.max)
        bar = '|' + '>'*round(self.progress*self.width) + ' ' * round(self.width*(1-self.progress)) + '|'
        pct = round(self.progress*100)
        str_pct = str(pct) if pct > 10 else ' '+str(pct)
        print(f'\r  {bar}   {str_pct}% : {key}', end='\r')
