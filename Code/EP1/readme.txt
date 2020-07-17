#
#	
#	+---------------------------------------+
#	|	Isadora Rodrigues Bisognin	|
#	|	Luísa Mendes Heise		|
#	+---------------------------------------+
#
# 
#    _____       _             _                        
#   / ____|     (_)           | |                       
#  | |  __ _   _ _  __ _    __| | ___   _   _ ___  ___  
#  | | |_ | | | | |/ _` |  / _` |/ _ \ | | | / __|/ _ \ 
#  | |__| | |_| | | (_| | | (_| |  __/ | |_| \__ \ (_) |
#   \_____|\__,_|_|\__,_|  \__,_|\___|  \__,_|___/\___/ 
#                                                      
#	Requisitos:
#		- matplotlib
#		- numpy
#
#	********************
#	* Método explícito *
#	********************
#
# Para fazer a simulação com o *método explícito*: 
# 	- digite "python3 ex1.py" no terminal.
# 	- escolha o item, o N e o lambda.
# 	- será mostrado o gráfico gerado
#	- será salvo o gráfico em uma pasta ex1/item onde item é "a1", "a2", "b" ou "c"
#	- dentro de um arquivo erro.txt (que estará na pasta ex1) constará o item, o lambda, o N, o M, o maior erro de truncamento e o maior erro absoluto para a simulação.
#	- O arquivo erro.txt conterá essas infos para todas as simulações feitas
#
#
#	**********************
#	* Métodos implícitos *
#	**********************
# Para fazer as simulações com os *métodos implícitos*
# 	- digite "python3 ex2.py" no terminal.
# 	- escolha o método (Euler ou Crank-Nicolson), item e o N.
# 	- será mostrado o gráfico gerado
#	- será salvo o gráfico em uma pasta ex2
#	- dentro de um arquivo erro.txt (que estará na pasta ex2) constará o item, o lambda, o N, o M, o maior erro de truncamento e o maior erro absoluto para a simulação.
#	- O arquivo erro.txt conterá essas infos para todas as simulações feitas
#
#
#  _____                                        _             /\/|       
# |  __ \                                      | |           |/\/        
# | |  | | ___   ___ _   _ _ __ ___   ___ _ __ | |_ __ _  ___ __ _  ___  
# | |  | |/ _ \ / __| | | | '_ ` _ \ / _ \ '_ \| __/ _` |/ __/ _` |/ _ \ 
# | |__| | (_) | (__| |_| | | | | | |  __/ | | | || (_| | (_| (_| | (_) |
# |_____/ \___/ \___|\__,_|_| |_| |_|\___|_| |_|\__\__,_|\___\__,_|\___/ 
#                                                         )_)                                                                                  
#
#
#		+-----------------+
# 		|TAREFA 1 - ex1.py|
#		+-----------------+
#
#-------------------------------------------------------------------------------------------
f(x, t, item, dx = None)

	Essa função é utilizada para cálculo da função f de um dado problema
	-----------
	Input:
		- x = posição na barra [0,1]
		- t = instante de tempo [0, T]
		- item = "a1", "a2", "b", "c" referente ao item da tarefa 1
		- dx = intervalo de discretização do comprimento da barra (usado no item c)
	-----------
	Output:
		- valor da função f
#-------------------------------------------------------------------------------------------
u_0(x, item)

	Essa função calcula a condição de contorno u_0 de um dado problema
	-----------
	Input:
		- x = posição na barra [0,1]
		- item = "a1", "a2", "b", "c" referente ao item da tarefa 1
	-----------
	Output:
		- valor da função u_0
#-------------------------------------------------------------------------------------------
g_1(t, item)

	Essa função calcula a condição de contorno g_1 de um dado problema
	-----------
	Input:
		- t = instante de tempo [0, T]
		- item = "a1", "a2", "b", "c" referente ao item da tarefa 1
	-----------
	Output:
		- valor da função g_1
#-------------------------------------------------------------------------------------------
g_2(t, item)

	Essa função calcula a condição de contorno g_2 de um dado problema
	-----------
	Input:
		- t = instante de tempo [0, T]
		- item = "a1", "a2", "b", "c" referente ao item da tarefa 1
	-----------
	Output:
		- valor da função g_2
#-------------------------------------------------------------------------------------------
real(i,k, dx, dt, item)

	Essa função calcula o valor da distribuição de temperatura na barra usando a expressão
	analítica. Trata-se, portanto, do valor "real"/verdadeiro de u(i,k)
	-----------
	Input:
		- i = indice da discretização do comprimento da barra [0,N]
		- k = indice da discretização do tempo [0,M]
		- dx = intervalo de discretização do comprimento da barra
		- dt = = intervalo de discretização do tempo
		- item = "a1", "a2" ou "b" referente ao item da tarefa 1
	-----------
	Output:
		- o valor analítico (verdadeiro) de u(i,k)
#-------------------------------------------------------------------------------------------
computa_us_v(solucao_atual, k, i, item, dx, dt)

	Essa função computa um u(i,k) para um dado k e i
	-----------
	Input:
		- solucao_atual = vetor solução referente ao intante de k-1
		- i = indice da discretização do comprimento da barra [0,N] 
		- k = indice da discretização do tempo [0,M]
		- dx = intervalo de discretização do comprimento da barra
		- dt = = intervalo de discretização do tempo
		- item = "a1", "a2", "b", "c" referente ao item da tarefa 1 
	-----------
	Output:
		- u(i,k) 
#-------------------------------------------------------------------------------------------
vectorize_us(N, M, dx, dt, item, lamb)

	Essa função realiza o loop principal de calculo de todos u(i,k) dentro da malha.
	Ela também calcula os erros absoluto e de truncamento.
	São feitos os plots da solução discretizada. 
	Além disso,é escrito um arquivo txt com o valor máximo do erro de truncamento e absoluto.

	Essa função não retorna nada. Apenas salva a figura com o gráfico e escreve o arquivo txt.
	-----------
	Input:
		- N = número de divisões do comprimento da barra
		- M = número de divisões do tempo
		- dx = intervalo de discretização do comprimento da barra
		- dt = = intervalo de discretização do tempo
		- item = "a1", "a2", "b", "c" referente ao item da tarefa 1 
		- lamb = fator lambda
		
#-------------------------------------------------------------------------------------------
main()
	Essa função recebe os inputs do usuário quanto ao item desejado, o N desejado e o
	lambda desejado e passa isso para a função vectorize_us
#-------------------------------------------------------------------------------------------
#
#		+-----------------+
# 		|TAREFA 2 - ex2.py|
#		+-----------------+
#
#------------------------------------------------------------------#
#					ITEM A			   #
#------------------------------------------------------------------#
l_d_l(A):
	
	Essa função calcula a decomposição LDL de uma matriz simétrica tridiagonal A
	-----------
	Input:
		- A = dois vetores representando a diagonal principal e a subdiagonal da
			matriz tridiagonal e simétrica em questão
	-----------
	Output:
		-   Dois vetores: com a diagonal da matriz D e subdiagonal das matrizes L e L^t
#-------------------------------------------------------------------------------------------
resolve_sis_L(L,b)
	
	Essa função resolve um sistema do tipo (m_L)x = b, em que m_L seria uma matriz com diagonal com 1s e com L como subdiagonal.
	
	-----------
	Input:
		- L: vetor representando a subdiagonal da matriz m_L
		- b: vetor representando o lado direito do sistema linear
	-----------
	Output:
		- vetor com a solução do sistema
#-------------------------------------------------------------------------------------------
resolve_sis_D(D,z)

	Essa função resolve um sistema do tipo (m_D)x = z, em que m_D seria uma matriz com diagonal com D como diagonal.
	
	-----------
	Input:
		- D: vetor representando a diagonal da matriz m_D
		- z: vetor representando o lado direito do sistema linear
	-----------
	Output:
		- vetor com a solução do sistema
#-------------------------------------------------------------------------------------------
resolve_sis_Lt(Lt,y)
	
	Essa função resolve um sistema do tipo (m_Lt)x = y, em que m_Lt seria uma matriz com diagonal com 1s e com Lt como subdiagonal.
	
	-----------
	Input:
		- Lt: vetor representando a subdiagonal da matriz m_Lt
		- y: vetor representando o lado direito do sistema linear
	-----------
	Output:
		- vetor com a solução do sistema
#-------------------------------------------------------------------------------------------
resolve_c_ldl(D,L,b)
	
	Essa função utiliza as funções resolve_sis_D, resolve_sis_L e resolve_sis_Lt para resolver um sistema tendo os vetores da decomposição LDLt de uma matriz tridiagonal simétrica e o lado direito da equação (vetor b)

	-----------
	Input:
		- D: vetor representando a diagonal da matriz D
		- L: vetor representando a subdiagonal da matriz L e Lt
		- b: vetor representando o lado direito do sistema linear
	-----------
	Output:
		- vetor com a solução do sistema
#------------------------------------------------------------------------------------------
resolucao(A,b)

	Essa função recebe uma matriz A e realiza sua decomposição LDL com uso da função l_d_l e depois resolve o sistema LDLt x = b com uso da função resolve_c_ldl.

	-----------
	Input:
		- A = dois vetores representando a diagonal principal e a subdiagonal da
			matriz tridiagonal e simétrica em questão
		- b: vetor representando o lado direito do sistema linear
	-----------
	Output:
		- vetor com a solução do sistema
#------------------------------------------------------------------#
#					ITENS B e C		   #
#------------------------------------------------------------------#
erro_de_truncamento(solucao_anterior, solucao_atual, dt, dx, k, ex, item)

	Essa função calcula o erro de truncamento para um dos métodos implícitos determinado para um dado instante de tempo
	
	------
	Input:
		- solucao_anterior: solução cujo instante seja calculado o erro de truncamento
		- solucao_atual: solução do passo seguinte a solução anterior
		- dt: intervalo de tempo
		- dx: intervalo de espaço da barra
		- k: k correspondente ao intervalo de tempo da solucao_anterior
		- ex: especificação do método implícito, crank-nicolson ou euler (0 ou 1)
		- item: especificação do item ("a1", "a2", "b" ou "c")

	--------
	Output:
		- vetor contendo os erros de truncamento
#-------------------------------------------------------------------------------------------
mtz_A(ex, lamb, N)

	Essa função monta a matriz do lado esquerdo do sistema linear para um dos métodos

	------
	Input:
		- ex: especificação do método implícito, crank-nicolson ou euler (0 ou 1)
		- lamb: fator lambda
		- N: número de divisões do comprimento da barra

	-------
	Output:
		- dois vetores, contendo a subdiagonal e a diagonal da matriz A formada

#-------------------------------------------------------------------------------------------
b_(ex, N, k, dt, dx, lamb, sol_anterior, item)

	Essa função calcula o vetor b relativa as condições de input

	------
	Input:
		- ex: especificação do método implícito, crank-nicolson ou euler (0 ou 1)
		- lamb: fator lambda
		- N: número de divisões do comprimento da barra
		- k: k correspondente ao intervalo de tempo
		- dt: intervalo de tempo
		- dx: intervalo de espaço da barra
		- solucao_anterior: solução cujo instante anterior
		- item: especificação do item ("a1", "a2", "b" ou "c")

	-------
	Output:
		- o vetor b (lado direito do sistema linear) correpondente

#-------------------------------------------------------------------------------------------
calcular_us_sistema(ex, item, N, M, lamb, dx, dt)

	Essa função realiza o loop principal de calculo de todos u(i,k) dentro da malha.
	Ela também calcula os erros absoluto e de truncamento.
	São feitos os plots da solução discretizada. 
	Além disso,é escrito um arquivo txt com o valor máximo do erro de truncamento e absoluto.

	Essa função não retorna nada. Apenas salva a figura com o gráfico e escreve o arquivo txt.
	-----------
	Input:
		- N = número de divisões do comprimento da barra
		- M = número de divisões do tempo
		- dx = intervalo de discretização do comprimento da barra
		- dt = = intervalo de discretização do tempo
		- item = "a1", "a2", "b", "c" referente ao item da tarefa 1 
		- lamb = fator lambda
		- ex: especificação do método implícito, crank-nicolson ou euler (0 ou 1)
	
#-------------------------------------------------------------------------------------------
main()

	Essa função recebe os inputs do usuário quanto ao item e método desejado, assim como o N e passa isso para a função calcular_us_sistema
