from motor import uncertainty

# Modelo de medição (cálculo de volume)
def model(L):
    L1 = L[0]; L2 = L[1]; L3 = L[2]
    return L1 * L2 * L3

# Dados de entrada
L1 = 11.58; L2 = 23.72; L3 = 33.16     # Medições
uL1 = 0.029; uL2 = 0.043; uL3 = 0.049  # Incertezas
rL1L2 = 0.000000000000000001; rL1L3 = 0.0000000000000001; rL2L3 = 0.00000000000000001 # Correlação
GL1 = 11; GL2 = 12; GL3 = 13 # Graus de liberdade

L = [L1,L2,L3]
u = [uL1,uL2,uL3]
r = [rL1L2, rL1L3, rL2L3]

incerteza = uncertainty(model,x=L,symbols=None,u=u,r=r)

uc = incerteza.combinated()
volume = model(L)

Gl_efetivo = (uc/volume)**4/(((uL1/L1)**4)/GL1+((uL2/L2)**4)/GL2+((uL3/L3)**4)/GL3)

print('medição = ' + str(volume))
print('Incerteza-padrão Combinada = ' + str(uc))
print('Gl = ' + str(Gl_efetivo))
print('K = ' + str(incerteza.expandida(PA=0.90, GL=Gl_efetivo)))
print('Incerteza Expandida = ' + str(incerteza.expandida(PA=0.90, GL=Gl_efetivo)*uc))

