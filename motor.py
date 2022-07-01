from casadi import MX,vertcat,jacobian, Function, dot
from math import sqrt
from scipy.stats import t

class uncertainty:

    def __init__(self, model, x, symbols, u, r):

        self.model = model      # modelo de medição
        self.u = u              # incertezas das grandezas de entrada
        self.symbols = symbols  # simbolos das grandezas de entrada
        self.r = r              # correlação entre as grandezas de entrada
        self.grandezas = x      # estimativa das grandezas de entrada

    def symbolic(self):
        
        # método para criar, no formato simbólico, as variáveis, o modelo de medição e a lei de propagação de incertezas (LPI)

        x = []
        for i in range(len(self.grandezas)):
            x = vertcat(x,MX.sym('q' + str(i)))

        # uA is always the first
        uA = MX.sym('uA')

        uB = []
        for i in range(len(self.u)-1):
            uB = vertcat(uB, MX.sym('uB'+ str(i+1)))

        u = vertcat(uA,uB)

        r = []
        for i in range(len(self.r)):
            r = vertcat(r, MX.sym('r' + str(i)))

        aux1 = 0; aux2 = 0; aux3 = 0

        self.measurementModel = self.model(x)

        for i in range(len(self.u)):
            aux1 = aux1 + dot(jacobian(self.measurementModel,x[i])**2,u[i]**2)

        k = 0
        for i in range(len(self.grandezas)-1):
            for j in range(i+1,len(self.grandezas),1):
                aux2 = aux2 + dot(dot(dot(dot(jacobian(self.measurementModel,x[i]), jacobian(self.measurementModel,x[j])),r[k]),u[i]),u[j])
                k = k + 1

        self.uCombinada_exec_independente = Function('Uc',[vertcat(x,u,r)],[aux1])     # Primeira parcela da LPI
        self.uCombinada_exec_dependente   = Function('Uc', [vertcat(x, u, r)], [aux2]) # Segunda parcela da LPI

    def resolution(self,R):

        # Método para calcular a incerteza associada à resolução

        return R/sqrt(12)

    def combinated(self):

        # Método para calcular a incerteza-padrão combinada
        self.symbolic()

        uc_dep = self.uCombinada_exec_dependente(vertcat(self.grandezas,self.u,self.r))   # Primeira parcela
        uc_ind = self.uCombinada_exec_independente(vertcat(self.grandezas,self.u,self.r)) # Segunda parcela

        uc = float((uc_ind + 2 * uc_dep) ** 0.5)

        return uc

    def expandida(self, PA, GL):

        # Método para calcular a incerteza expandida

        k = -t.ppf((1-PA)/2,GL)
        uc = self.combinated()
        U = uc*k

        return k



