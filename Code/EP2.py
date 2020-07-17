#Começamos importando as funções do EP1 assim como o numpy e re (para limpar o txt)
import numpy as np
from EP1 import ex2
import re
import warnings
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")

#=========================================================
#  Funções Principais
#

def main():
    """
    Essa é a função main do EP, aqui vamos pedir os inputs do usuário e direcionar esses inputs para a função de simulação
    """
    item = input("Qual item deseja fazer a simulação (a, b, c ou d)? ")
    assert item in ["a","b", "c", "d"]

    if item == "a":
        pks = [0.35]
        #N = 128
        N = int(input("Digite o N: "))
        h = 1/N
        simulacao(pks, N, h, "a")

    if item == "b":
        pks = [0.15,0.3,0.7,0.8]
        #N = 128
        N = int(input("Digite o N: "))
        h = 1/N
        simulacao(pks, N, h, "b")
    
    if item == "c":
        pks, uT = leitura_txt()
        #Ns = [128, 256, 512, 1024, 2048]
        N = int(input("Digite o N: "))
        step = int(2048/N)
        uT_step = [uT[i] for i in range(0,2049,step)][1:N]
        h = 1/N
        print("Por favor aguarde...")
        simulacao(pks, N, h, "c", uT_step)
        grafico_uT(uT, uT_step, N)

    
    if item == "d":
        pks, uT = leitura_txt()
        N = int(input("Digite o N: "))
        step = int(2048/N)
        uT_step = [uT[i] for i in range(0,2049,step)][1:N]
        for i in range(len(uT_step)):
        #multiplicação por número aleatório
            multiplicador = 1 + 0.01*np.random.uniform(low=-1, high=1)
            uT_step[i] = multiplicador*uT_step[i]

        #Ns = [128, 256, 512, 1024, 2048]
        h = 1/N
        print("Por favor aguarde...")
        simulacao(pks, N, h, "d", uT_step)
        grafico_uT(uT, uT_step, N)


def simulacao(pks, N, h, item, uT=None):
    """
    Essa função faz a rotina de simulação:
        - Determinação dos valores de uk(T, x_i) tendo os forçantes pontuais, condições de contorno nulas e utilizando o método de crank-nicolson
        - Determinação dos valores de uT(x_i), quando estes são função de uk(T, x_i)
        - Determinação dos valores de a por meio do método dos mínimos quadrados
    Ao final, são impressos os valores de a e o erro E2, quando for um caso de item c ou d.

    Keyword arguments:
        - pks: lista contendo os valores dos forçantes pontuais
        - N: número referente ao tamanho do intervalo dx = 1/N da discretização da barra
        - h: argumento referente à zona de atuação do forçante pontual, é um argumento da função g
        - item: refere-se ao item do enunciado, ou seja, as condições da simulação
        - uT: vetor contendo uT no caso em que ele é pré-determinado
    """
#calculo dos uk[T]
    uks = vetores_uk(pks, N, h)

#nos casos em que o uT é uma função dos uks, definimos o vetor uT
    if item == "a":
        uT = 7*uks[0]
    if item == "b":
        uT = 2.3*uks[0] + 3.7*uks[1] + 0.3*uks[2] + 4.2*uks[3]
#calculamos o vetor com os coeficientes a com o mmq
    a_vector = resolve_mmq(uks, uT)
#Imprimindo as respostas
    print(f"N = {N}")
    for i,a in enumerate(a_vector):
        print(f"a({i}) = {a}")
#Erro quadrático, quando couber (itens c e d)
    if item in ["c", "d"]:
        E2 = e2(1/N, a_vector, uks, uT, N)
        print(f"O erro E2 = {E2}")

def grafico_uT(uT_real, uT, N):
    """
    Essa função plota as soluções real e utilizada para efeitos de comparação

    Keyword arguments:
        - uT_real: vetor contendo uT verdadeiro 
        - uT: vetor contendo uT utilizado 
        - N: número referente ao tamanho do intervalo dx = 1/N da discretização da barra
    """
    plt.figure(figsize=(10,5))
    #condição de contorno
    dx = 1/N
    x_barra = np.array([dx*n for n in range(N+1)])
    x_barra_real = np.array([(1/2048)*n for n in range(2049)])
    uT = np.concatenate([[0], uT, [0]])
    plt.plot(x_barra_real, uT_real, label="uT real", color='black')
    plt.plot(x_barra, uT, label="uT simulação", color='red', linestyle="dashed")
    plt.xlabel("Comprimento da barra")
    plt.ylabel("Distribuição de temperatura em T")
    plt.title(f"u_T(x) verdadeiro x uT(x) simulação com N={N}")
    plt.legend()
    plt.show()

#=========================================================
#leitura do txt
def leitura_txt():
    """
    Essa função faz a leitura do arquivo txt contendo os forçantes pontuais e uT no arquivo teste.txt 
    """

    f = open("teste.txt", "r")
    #a primeira linha contem as fontes
    fontes_raw = f.readline()
    #aplicamos uma regex para tirar os espaços e enteres, habilitando a string a ser convertida em um float
    #lendo o vetor com os forçantes pontuais
    pks = [float(re.sub("[^0-9,.,E,-]", "", fonte)) for fonte in fontes_raw.split("       ")]
    uT = []
    #lendo as linhas
    for i in range(2049):
        ut = f.readline()
        ut = float(re.sub("[^0-9,.,E,-]", "", ut))
        uT.append(ut)

    return pks, uT

def e2(dx, a, uks, uT, N):
    """ 
    Essa função calcula o erro E2

    Keyword arguments:
        - dx: intervalo de discretização do comprimento da barra
        - a: vetor contendo os valores de a calculados
        - uks: vetores uks calculados
        - uT: vetor contendo uT
        - N: número referente ao tamanho do intervalo dx = 1/N da discretização da barra
    
    Output:
        - E2: erro E2 
    """
    nf = len(uks)
    soma_quadrados = 0
    for i in range(N-1):
        soma_ak_uk = sum([a[k]*uks[k][i] for k in range(nf)])
        soma_quadrados += (uT[i] - soma_ak_uk)**2
    E2 = (dx*soma_quadrados)**(1/2)
    return E2

#==========================================================
#As funções para a geração da função f e a própria função f
def g_kh(x, pk, h):
    """
    Cálculo da função g(x), utilizada no cálculo de f(t,x)

    Keyword arguments:
        - x: coordenada na barra
        - pk: forçante pontual
        - h: intervalo de atuação do forçante pontual

    Output:
        - g: valor da função g(x)
    """
    if (x >= pk - (h/2)) and (x <= pk + (h/2)):
        g = (1/h)
    else:
        g = 0
    return g

def r(t):
    """
    Cálculo da função r(t), utilizada no cálculo de f(t,x)

    Keyword arguments:
        - t: instante de tempo [0,T]

    Output:
        - r: valor da função r(t)
    """
    return 10*(1 + np.cos(5*t))

def f_tx(t, x, pk, h):
    """
    Cálculo da função f(t,x) para os diversos calculos de crank nicolson

    Keyword arguments:
        - t: instante de tempo [0,T]
        - x: coordenada na barra [0,1]
        - pk: forçante pontual
        - h: intervalo de atuação do forçante pontual

    Output:
        - f: valor da função f(t,x)
    """
    r_t = r(t)

    g_x = g_kh(x, pk, h)
    f = r_t*g_x
    return f
#============================================================

def b_(k, pk, h, sol_anterior, lamb, dt, dx, N):
    """
    Vetor b do lado direiro da equação do método de crank nicolson, já considerando que as condições de 
    contorno g1, g2 e u0 são nulas

    Keyword arguments:
        - k: índice do instante de tempo [0,M]
        - pk: forçante pontual
        - h: intervalo de atuação do forçante pontual
        - sol_anterior: solução referente ao instante de tempo de índice k-1
        - lamb: parâmetro lambda 1/dx
        - dt: intervalo de discretização do tempo T/M
        - dx: intervalo de discretização da barra 1/N
        - N: número de discretizações da barra (N=M)

    Output:
        - b: valor do lado direito da equação do método de crank nicolson para o instante de índice k
    """
    b = np.zeros(N-1)
    t = dt*k
    b[0] = sol_anterior[0] + (lamb/2)*(-2*sol_anterior[0] + sol_anterior[1]) + (dt/2)*(f_tx((k-1)*dt, dx, pk, h) + f_tx(t, dx, pk, h))
    for i in range(1,N-2):
        x=(i+1)*dx
        b[i] = sol_anterior[i] + (lamb/2)*(sol_anterior[i-1] - 2*sol_anterior[i] + sol_anterior[i+1]) + (dt/2)*(f_tx(t, x, pk, h) + f_tx(dt*(k-1), x, pk, h))
    b[N-2] = sol_anterior[N-2] + (lamb/2)*(sol_anterior[N-3] - 2*sol_anterior[N-2]) + (dt/2)*(f_tx(dt*(k-1), dx*(N-1), pk, h) + f_tx(t, dx*(N-1), pk, h))
    return b  

#solve crank nicolson

def crank_nicolson(A, M, N, lamb, dt, dx, h, pk):
    """
    Calculo de um dos vetores de uk(T) referente a um dos forçantes pontuais (pk).
    Para isso, são calculados todos os vetores uk(t) para os instantes de tempo anteriores a T.
    Isso porque o crank nicolson requer o vetor anterior para calcular o atual

    Keyword arguments:
        - A: matriz do lado do crank-nicolson
        - N: número de discretizações da barra (N=M)
        - M: número de discretizações do tempo 
        - lamb: parâmetro lambda 1/dx
        - dt: intervalo de discretização do tempo T/M
        - dx: intervalo de discretização da barra 1/N
        - h: intervalo de atuação do forçante pontual
        - pk: forçante pontual

    Output:
        - u: valor da distribuição de temperatura com f dependendo do forçante pontual uk no instante t=T (exclui as posições 0 e N!)
    """
    #A primeira solução determinada é u0, que eh uma condição de controno nula:
    sol_anterior = np.zeros(N+1)
    for k in range(1,M+1):
        b = b_(k, pk, h, sol_anterior, lamb, dt, dx, N)
        u = ex2.resolucao(A,b)
        sol_anterior = u
    u = np.array(u).reshape(-1)
    #Depois do calculo do ultimo U
    return u

def vetores_uk(pks, N, h):
    """
    Calculo dos vetores uks(T) referente aos forçantes pontuais (pks).
    Para cada forçante pontual, é calulado o vetor uk(T) correspondente via crank nicolson

    Keyword arguments:
        - pks: vetor com forçantes pontuais
        - N: número de discretizações da barra (N=M)
        - h: intervalo de atuação do forçante pontual

    Output:
        - uks: lista com uks: valor da distribuição de temperatura com f dependendo do forçante pontual 
        uk no instante t=T para cada um dos pks 
    """
    T = 1
    M = N
    dt, dx = T/M, 1/N
    lamb = N
    #incializando a lista que conterá todos os vetores uk
    uks = []
    #crank nicolson para obter os uk(T,x), chamados de uk
    A = ex2.mtz_A(1,lamb,N)
    for pk in pks:
        uk = crank_nicolson(A, M, N, lamb, dt, dx, h, pk)
        uks.append(uk)        
    return uks
#===========================================================
# Sistema normal MMQ, LDLt e eliminação de gauss

def mtz_normal(uks):
    """
    Calculo da matriz normal do MMQ

    Keyword arguments:
        - uks: lista com uks: valor da distribuição de temperatura com f dependendo do forçante pontual 
        uk no instante t=T para cada um dos pks 

    Output:
        - mtz: matriz normal do MMQ
    """
    #a dimensao da matriz normal é nfXnf 
    #ou seja len(uks)Xlen(uks)
    mtz = []
    uks = np.array(uks)
    for i in range(len(uks)):
        #montando as linhas
        linha = []
        for j in range(len(uks)):
            #o produto interno tal como definido é um vetor vezes o outro transposto
            elemento = uks[i].T@uks[j]
            linha.append(elemento)
        mtz.append(linha)
    return mtz

#decomposição LDLt
def LDLt(matriz):
    """
    Calculo da decomposição LDLt da matriz (será inoutada a matriz normal, mas essa função é genérica)

    Keyword arguments:
        - matriz: matriz para a qual deseja-se fazer a decomposição

    Output:
        - L: matriz triangular inferior referente a decomposição LDLt
        - D: vetor contendo os elementos da diagonal da matriz diagonal da decomposição LDLt
    """
    n = len(matriz)
    D = np.zeros(n)
    L = np.zeros((n,n))
    for i in range(n):
        v = [L[i][j]*D[j] for j in range(i)]
        soma = sum([L[i][j]*v[j] for j in range(i)])
        D[i] = matriz[i][i] - soma
        for j in range(i+1, n):
          soma2 = 0
          for k in range(i):
            soma2 += L[j][k]*v[k]
          L[j][i] = (matriz[j][i] - soma2)/D[i] 
    #diagonal com 1s
    for i in range(n):
        L[i][i] = 1
    return L, D


def eliminacao_de_gauss(matriz, resposta):
    """
    Calculo de x em: matriz@x = resposta, mais comumente representado como Ax=b
    Esse processo é feito via eliminação de gauss

    Keyword arguments:
        - matriz: matriz do lado esquerdo do sistema linear
        - reposta: vetor correspondente ao lado direito do sistema linear

    Output:
        - x: resposta do calculo do sistema linear
    """
    #Essa função pressuopõe que os inputs sejam listas e não np.arrays!
    if not(isinstance(matriz, list)):
        matriz = matriz.tolist()

    if not(isinstance(resposta, list)):
        resposta = resposta.tolist()

    n = len(matriz)
    m = len(matriz[0]) #número de colunas
    assert n == m, "Esse sistema linear não tem nenhuma solução ou infinitas ou esta superespecificado"
#se o número de linhas e de colunas for diferente não faremos o processo
#adicionando a resposta na matriz 
    i = 0
    for linha in matriz:
        linha.append(resposta[i])
        i += 1
#pivotamento
    for p in range(n):
        j = p
        solucionavel = False
        while j<n and (not solucionavel):
            if matriz[p][j] != 0:
                linha = p
                solucionavel = True
        else:
            pass
        if linha != p:
            copia_matriz = matriz[:]
            matriz[p] = copia_matriz[linha]
            matriz[linha] = copia_matriz[p]

#resolvendo o sistema, fazendo a eliminiação de gauss em si
        for j in range(p+1,n):
            coeficiente = float(matriz[j][p]) / matriz[p][p]
            for m in range(p, n+1):
                matriz[j][m] -=  coeficiente * matriz[p][m]

#preenchendo provisoriamente o velor solucação com zeros
    solucao = [0]*n
    assert matriz[n-1][n-1] != 0, "Não exite uma solução única para o sistema" 
    solucao[n-1] =float(matriz[n-1][n])/matriz[n-1][n-1]
    for i in range (n-1,-1,-1):
        soma = 0
        for j in range(i+1,n):
            soma = soma  + float(matriz[i][j])*solucao[j]
        solucao[i] = float(matriz[i][n] - soma)/matriz[i][i]
    return solucao


def resolve_mmq(uks, uT):
    """
    Essa função calcula os vetores a via método dos mímimos quadrados, para isso são seguidos os seguintes passos:
        1. Calculo da matriz normal (com os uks)
        2. Cálculo do vetor do lado direito do sistema normal (b_normal)
        3. Cálculo da decomposição LDLt da matriz normal
        4. Cálculos de sistemas intermediários até obter a.

    Keyword arguments:
        - uks: lista com uks: valor da distribuição de temperatura com f dependendo do forçante pontual 
        uk no instante t=T para cada um dos pks 
        - uT:  vetor contendo uT 

    Output:
        - a: vetor contendo os valores dos coeficientes das fontes
    """

    normal = mtz_normal(uks)
    ##
    if not(type(uT) is np.ndarray):
        uT = np.array(uT)

    if not(type(uks) is np.ndarray):
        uks = np.array(uks)
    ##
    b_normal = [np.array([uT.T@uk]).reshape(-1)[0] for uk in uks]
    if len(uks) == 1:
        a = [b_normal[0]/normal[0][0]]
    else:
        L, D = LDLt(normal)
        #Temos Lw = b => Dy = w => Lta = y
        w = eliminacao_de_gauss(L, b_normal)
        y = ex2.resolve_sis_D(D,w).reshape(-1)
        a = eliminacao_de_gauss(L.T, y)
    return a

####################
main()