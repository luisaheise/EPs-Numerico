import numpy as np
import matplotlib.pyplot as plt
import ex1
import os
#----------------------------------------------------------------------------------
# Pastas
if not os.path.exists('ex2'):
    os.makedirs('ex2')
#----------------------------------------------------------------------------------
#-------------- 2 A ---------------------------------------------------------------
#----------------------------------------------------------------------------------
def l_d_l(A):
    diagonal, subdiagonal = A[0], A[1]

    #incializando os vetores da decomposição LDL
    D = np.zeros((len(diagonal),1))
    L = np.zeros((len(diagonal)-1,1))

    #incializando
    d_a = diagonal[0]
    D[0] = d_a
    l_a = subdiagonal[0]/d_a
    L[0] = l_a

    for i in range(1,len(diagonal)):
        if i == (len(diagonal)-1):
            d_d = diagonal[i] - d_a*(l_a**2)
            D[i]= d_d
        else:
            d_d = diagonal[i] - d_a*(l_a**2)
            l_d = subdiagonal[i]/d_d
            D[i], L[i] = d_d, l_d
            d_a, l_a = d_d, l_d
    
    return D, L
#----------------------------------------------------------------------------------
def resolve_sis_L(L,b):
    z = np.zeros((len(L)+1,1))
    z[0] = b[0]
    for i in range(1,len(z)):
        z[i] = b[i] - L[i-1]*z[i-1]
    return z
#----------------------------------------------------------------------------------
def resolve_sis_D(D,z):
    y = np.zeros((len(D),1))
    for i in range(len(y)):
        y[i] = z[i]/D[i]
    return y
#----------------------------------------------------------------------------------
def resolve_sis_Lt(Lt,y):
    x = np.zeros((len(Lt)+1,1))
    x[len(x)-1] = y[len(x)-1]
    for i in range(len(x)-2,-1,-1):
        x[i] = y[i] - Lt[i]*x[i+1]
    return x
#----------------------------------------------------------------------------------
def resolve_c_ldl(D,L,b):
    z = resolve_sis_L(L,b)
    y = resolve_sis_D(D,z)
    x = resolve_sis_Lt(L,y)
    return x
#----------------------------------------------------------------------------------
def resolucao(A,b):
    D, L = l_d_l(A)
    u = resolve_c_ldl(D,L,b)
    return u
#----------------------------------------------------------------------------------
#-------------- 2 B e C -----------------------------------------------------------
#----------------------------------------------------------------------------------
def erro_de_truncamento(solucao_anterior, solucao_atual, dt, dx, k, ex, item):
    assert solucao_anterior.shape == solucao_atual.shape, "Erro de dimensão"
    vetor_erro_de_truncamento =[]
    if ex == 0:
        for i in range(len(solucao_anterior)):
            if i == 0:
                erro_trunk_i = 0
            elif i == (len(solucao_anterior) - 1):
                erro_trunk_i = 0
            else:
                p1 = (solucao_atual[i] - solucao_anterior[i])/dt
                p2 = (solucao_atual[i-1] - 2*solucao_atual[i] + solucao_atual[i+1])/(dx**2)
                erro_trunk_i = p1 - p2 - ex1.f(dx*i, (k+1)*dt, item, dx)
            vetor_erro_de_truncamento.append(erro_trunk_i)
    if ex == 1:
        for i in range(len(solucao_anterior)):
            if i == 0:
                erro_trunk_i = 0
            elif i == (len(solucao_anterior) - 1):
                erro_trunk_i = 0
            else:
                p1 = solucao_atual[i] - solucao_anterior[i]
                p2_1 = (solucao_atual[i-1] - 2*solucao_atual[i] + solucao_atual[i+1])
                p2_2 = (solucao_anterior[i-1] - 2*solucao_anterior[i] + solucao_anterior[i+1])
                p2 = p2_1 + p2_2
                p3 = ex1.f(dx*i, k*dt, item, dx) + ex1.f(dx*i, (k+1)*dt, item, dx)
                erro_trunk_i = (p1/dt) - (p2/(2*(dx**2))) - p3/2
            vetor_erro_de_truncamento.append(erro_trunk_i)
    return np.array(vetor_erro_de_truncamento)
#----------------------------------------------------------------------------------
def mtz_A(ex, lamb, N):
    if ex == 0:
        diagonal = np.repeat((1+2*lamb), N-1)
        subdiagonal = np.repeat((-lamb), N-2)
        A = np.array((diagonal,subdiagonal))
    if ex == 1:
        diagonal = np.repeat((1+lamb), N-1)
        subdiagonal = np.repeat((-lamb/2), N-2)
        A = np.array((diagonal,subdiagonal))
    return A
#----------------------------------------------------------------------------------
def b_(ex, N, k, dt, dx, lamb, sol_anterior, item):
    b = np.zeros(N-1)
    t = dt*k
    if ex == 0:
        b[0] = sol_anterior[0] + dt*ex1.f(dx,t,item,dx) + lamb*ex1.g_1(t,item)
        for i in range(1,N-2):
            x = (i+1)*dx
            b[i] = sol_anterior[i] + dt*ex1.f(x,t,item,dx)
        b[N-2] = sol_anterior[N-2] + dt*ex1.f((N-1)*dx,t,item,dx) + lamb*ex1.g_2(t,item)
    if ex == 1:
        b[0] = sol_anterior[0] + (lamb/2)*(ex1.g_1(dt*(k-1),item) - 2*sol_anterior[0] + sol_anterior[1]) + (dt/2)*(ex1.f(dx,(k-1)*dt,item,dx) + ex1.f(dx,t,item,dx)) + (lamb/2)*ex1.g_1(t,item)
        for i in range(1,N-2):
            x=(i+1)*dx
            b[i] = sol_anterior[i] + (lamb/2)*(sol_anterior[i-1] - 2*sol_anterior[i] + sol_anterior[i+1]) + (dt/2)*(ex1.f(x,(k-1)*dt,item,dx) + ex1.f(x,t,item,dx))
        b[N-2] = sol_anterior[N-2] + (lamb/2)*(sol_anterior[N-3] - 2*sol_anterior[N-2] + ex1.g_2((k-1)*dt,item)) + (dt/2)*(ex1.f(dx*(N-1),(k-1)*dt,item,dx) + ex1.f(dx*(N-1),t,item,dx)) + (lamb/2)*ex1.g_2(t,item)
    return b
#----------------------------------------------------------------------------------
def calcular_us_sistema(ex, item, N, M, lamb, dx, dt):

    v_u_0 = np.array([ex1.u_0(dx*i,item) for i in range(N+1)])
    #plotando
    if item =="a1"or item=="a2" or item == "b":
        max_erro_absoluto = 0
        max_erro_truncamento = 0
        figure, axes = plt.subplots(nrows=1, figsize=(10,10))
        axes.plot(dx*np.array(range(N+1)),v_u_0, 'o-', label = 't=0.0') 
    if item =="c":
        max_erro_absoluto = "NA"
        max_erro_truncamento = 0
        figure, axes = plt.subplots(nrows=1, figsize=(10,10)) 
        axes.plot(dx*np.array(range(N+1)),v_u_0, 'o-', label = 't=0.0')

    A = mtz_A(ex, lamb, N)

    #primeiro erro
    if item =="a1" or item =="a2" or item == "b":
        real__ = [ex1.real(i,0, dx, dt, item) for i in range(N+1)]
        erro = real__ - v_u_0
        if np.max(np.absolute(erro)) > max_erro_absoluto:
            max_erro_absoluto = np.max(np.absolute(erro))

    sol_anterior = v_u_0[1:N]
    n = 1

    #incia o loop
    for k in range(1,M+1):
        t_k = dt*k
        b = b_(ex, N,k,dt,dx,lamb,sol_anterior, item)
        u = resolucao(A,b)

        #incluindo condições de contorno na solução atual
        solution_plot = np.concatenate(([ex1.g_1(dt*k,item)],np.reshape(u,(N-1))))
        solution_plot = np.concatenate((solution_plot, [ex1.g_2(dt*k,item)]))

        #erro absoluto para as equações com solução analítica
        if item =="a1" or item =="a2" or item == "b":
            real__ = [ex1.real(i,k,dx,dt,item) for i in range(N+1)]
            erro = real__ - solution_plot
            if np.max(np.absolute(erro)) > max_erro_absoluto:
                max_erro_absoluto = np.max(np.absolute(erro))

        #incluindo condições de contorno na solução anterior
        solution_ant = np.concatenate(([ex1.g_1(dt*(k-1),item)],np.reshape(sol_anterior,(N-1))))
        solution_ant = np.concatenate((solution_ant, [ex1.g_2(dt*(k-1),item)]))

        #erro de trucamento da solução anterior
        erro_de_truncamento_k_anterior = erro_de_truncamento(solution_ant, solution_plot, dt, dx, k-1, ex, item)
        if np.max(np.absolute(erro_de_truncamento_k_anterior)) > max_erro_truncamento:
                max_erro_truncamento = np.max(np.absolute(erro_de_truncamento_k_anterior))
        
        #plotando de 0.1 em 0.1
        if (((t_k >= (n*0.1) - (0.0000000001)) and (t_k <= (n*0.1) + (0.0000000001))) and n<10) or k == M:   
            n+=1
            axes.plot(dx*np.array(range(N+1)), solution_plot, 'o-', label = ('t='+str(round(dt*k,1))))
        sol_anterior = u

    axes.set_title("Iterações da discretização", y=1.03, fontsize=20)
    axes.set_xlabel("Comprimento da barra")
    axes.set_ylabel("Distribuição de calor")
    axes.legend()
    if ex == 0:
        exname = "euler"
    if ex == 1:
        exname = "crank_nicolson"
    path_fig = "ex2/" + str(exname)+"_"+str(item)+"_N="+str(N)+".png"
    plt.savefig(path_fig)

    f_path = f"ex2/erro_{exname}.txt"
    registro = open(f_path, "a")
    content = f"{item},{N},{M},{lamb},{max_erro_absoluto},{max_erro_truncamento}\n"
    registro.write(content)
    registro.close()

    plt.show()
#----------------------------------------------------------------------------------
def main():
    T = 1
    ex = int(input("Você quer simular o método de Euler (digite 0) ou Crank-Nicolson (digite 1)? "))
    #assert ex == 0 or ex == 1, "O método escolhido deve ser 0 ou 1"
    item = str(input("formulas do ex 1a1 (digite a1) ou 1a2 (digite a2) ou 1b(digite b) ou 1c(digite c)? "))
    assert(item == "a1" or item == "a2" or item == "b" or item == "c"), "escolha um item válido"
    N = int(input("Digite o número de subdivisões do espaço: "))
    dx = float(1/N)
    dt = dx
    lamb = 1/dx
    M = int(T/dt)
    calcular_us_sistema(ex, item, N, M, lamb, dx, dt)         
#----------------------------------------------------------------------------------
if __name__ == "__main__":
    main()