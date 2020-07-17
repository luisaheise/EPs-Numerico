import numpy as np
import matplotlib.pyplot as plt
import os
#----------------------------------------------------------------------------------
# Pastas
directories = ['ex1', 'ex1/a1', 'ex1/a1', 'ex1/b', 'ex1/c']
for directory in directories:
    if not os.path.exists(directory):
        os.makedirs(directory)
#----------------------------------------------------------------------------------
def f(x, t, item, dx = None):
    if item == "a1":
        return ((10*(x**2)*(x-1)) - (60*x*t) +(20*t))
    if item == "a2":
        p1 = 10*np.cos(10*t)*(x**2)*((1-x)**2)
        p2 = (1 + np.sin(10*t))*(12*(x**2)-12*x+2)
        return (p1 - p2)
    if item == "b":
        expo = np.exp(t-x)
        cosseno = np.cos(5*t*x)
        seno = np.sin(5*t*x)
        p1 = 25*expo*(t**2)*cosseno
        p2 = 10*expo*t*seno
        p3 = 5*expo*x*seno
        return (p1 - p2 - p3)
    if item == "c":
        p = 0.25
        h = dx
        h__ = h/2
        r_t = 10000*(1-2*(t**2))
        if p-h__ <= x <= p+h__:
            g_hx = 1/h
        else:
            g_hx = 0
        return r_t*g_hx
#-----------------------------------------------------------------------------------
def u_0(x, item):
    if item == "a1" or item == "c":
        return 0
    if  item == "a2":
        return (x**2)*((1-x)**2)
    if item == "b":
        return np.exp(-x)
#-----------------------------------------------------------------------------------
def g_1(t, item):
    if item == "a1" or item == "a2" or item == "c":
        return 0
    if item == "b":
        return np.exp(t)
#-----------------------------------------------------------------------------------
def g_2(t, item):
    if item == "a1" or item == "a2" or item == "c":
        return 0
    if item == "b":
        return np.exp(t-1)*np.cos(5*1*t)
#-----------------------------------------------------------------------------------
def real(i,k, dx, dt, item):
    x_i = dx*i
    t_k = dt*k
    if item == "a1":
        return 10*(t_k)*(x_i**2)*(x_i-1)
    if item == "a2":
        return (1+np.sin(10*t_k))*(x_i**2)*((1-x_i)**2)
    if item == "b":
        return np.exp(t_k-x_i)*np.cos(5*t_k*x_i)
#----------------------------------------------------------------------------------
def erro_de_truncamento(anterior_atual, k, item, dx, dt, N):
    erro_de_truncamento = []
    for i in range(anterior_atual[0].shape[0]):
        if i == 0:
            err_u = 0
        elif i == N:
            err_u = 0
        else:
            err_u = ((anterior_atual[1][i] - anterior_atual[0][i])/dt) - ((anterior_atual[0][i-1] - 2*anterior_atual[0][i] + anterior_atual[0][i+1])/(dx**2)) - f(i*dx,dt*(k-1),item,dx)
        erro_de_truncamento.append(err_u)
    return erro_de_truncamento
#----------------------------------------------------------------------------------
def computa_us_v(solucao_atual, k, i, item, dx, dt):
    x_i = dx*i
    t_k_1 = dt*(k-1)
    fracao = (solucao_atual[i-1] - 2*solucao_atual[i] + solucao_atual[i+1])/(dx**2)
    sol = solucao_atual[i] + dt*(fracao + f(x_i, t_k_1, item, dx))
    return sol
#----------------------------------------------------------------------------------
def vectorize_us(N, M, dx, dt, item, lamb):
    anterior_atual = np.zeros((2, N+1))
    #u_0
    anterior_atual[1] = [u_0(dx*i, item) for i in range(N+1)]
    #plot linha de u_0
    if item =="a1"or item=="a2" or item == "b":
        max_erro_absoluto = 0
        max_erro_truncamento = 0
        figure, axes = plt.subplots(nrows=1, figsize=(10,10))
        axes.plot(dx*np.array(range(N+1)), anterior_atual[1], 'o-', label = 't=0.0') 
    if item =="c":
        max_erro_truncamento = 0
        max_erro_absoluto = "NA"
        figure, axes = plt.subplots(nrows=1, figsize=(10,10))
        axes.plot(dx*np.array(range(N+1)), anterior_atual[1], 'o-', label = 't=0.0')
    #primeiro plot de erro
    if item =="a1" or item =="a2" or item == "b":
            real__ = [real(i,0,dx,dt,item) for i in range(N+1)]
            erro = real__ - anterior_atual[1]
            if np.max(np.absolute(erro)) > max_erro_absoluto:
                max_erro_absoluto = np.max(np.absolute(erro))

    n = 1
    #calculos de u
    for k in range(1,M+1):
        t_k = dt*k
        anterior_atual[0] = anterior_atual[1]
        #g_1
        anterior_atual[1][0] = g_1(t_k, item)
        #g_2
        anterior_atual[1][N] = g_2(t_k, item)
        #calculo
        anterior_atual[1][1:N] = np.array([computa_us_v(anterior_atual[0], int(k), int(i), item, dx, dt) for i in range(1,N)])
        #erro absoluto para as equações com solução analítica
        if item =="a1" or item =="a2" or item == "b":
            real__ = [real(i,k,dx,dt,item) for i in range(N+1)]
            erro = real__ - anterior_atual[1]
            if np.max(np.absolute(erro)) > max_erro_absoluto:
                max_erro_absoluto = np.max(np.absolute(erro))
        #erro de trucamento da solução anterior
        erro_de_truncamento_k_anterior =  erro_de_truncamento(anterior_atual, k, item, dx, dt, N)
        if np.max(np.absolute(erro_de_truncamento_k_anterior)) > max_erro_truncamento:
                max_erro_truncamento = np.max(np.absolute(erro_de_truncamento_k_anterior))
        #plotando em casos de multiplos de 0.1
        if (((t_k >= (n*0.1) - (0.000000001)) and (t_k <= (n*0.1) + (0.000000001)))and n<10 ) or (k == M):
            axes.plot(dx*np.array(range(N+1)), anterior_atual[1], 'o-', label = ('t='+str(round(t_k,1))))
            n+=1

    axes.set_title("Iterações da discretização", y=1.03, fontsize=20)
    axes.set_xlabel("Comprimento da barra")
    axes.set_ylabel("Distribuição de calor")
    axes.legend()

    path_fig = f"ex1/{item}/N_{N}_lamb_{1000*lamb}.png"
    plt.savefig(path_fig)

    f_path = "ex1/erro.txt"
    registro = open(f_path, "a")
    content = f"{item},{N},{M},{lamb},{max_erro_absoluto},{max_erro_truncamento}\n"
    registro.write(content)
    registro.close()

    plt.show()
#-----------------------------------------------------------------------------------
def main():
    T = 1
    item = input("Digite a1 para o primeiro subitem do item A, a2 para o segundo subitem do item A, b para o item B e c para o item C: ")
    assert item == "a1" or item == "a2" or item =="b" or item == "c", "Item invalido!"
    N = int(input("Digite o número de subdivisões do espaço (N): "))
    dx = 1/N
    lamb = float(input("Digite o λ: "))
    dt = lamb*(dx**2)
    M = int(T/dt)
    vectorize_us(N, M, dx, dt, item, lamb)

if __name__ == "__main__":
    main()