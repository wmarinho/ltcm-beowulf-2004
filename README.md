SIMULAÇÃO NUMÉRICA DA CONDUÇÃO DE CALOR TRIDIMENSIONAL EM UM CLUSTER BEOWULF
===

Resumo: O presente trabalho apresenta os recursos envolvidos na implementação de um cluster classe Beowulf para resolver problemas de transferência de calor tridimensional. O fenômeno simulado refere-se à condução de calor no interior de um sólido homogêneo e isotrópico, onde o domínio computacional - discretizado por diferenças finitas - pode ser dividido entre vários processadores em um ambiente distribuído. Além disso, são analisados tópicos importantes em processamento paralelo como características de hardware e software, decomposição do domínio, balanceamento de carga, problemas nas regiões de interface entre os subdomínios, troca de informações nas fronteiras, speedup e eficiência. Os resultados de speedup e eficiência são explorados para avaliar a performance do programa paralelo quando comparado a um programa serial. A curva de speedup é obtida a partir de um problema de tamanho fixo onde o número de processadores foi sendo variado, medindo-se o tempo para cada configuração. Além disso, a sua análise auxilia na determinação do número máximo de processadores possível de ser utilizado em um determinado programa paralelizado, otimizando o uso dos recursos computacionais. O trabalho mostra, também, o comportamento da curva de speedup em função do refinamento da malha e a importante relação da granularidade com a eficiência do software. O código computacional foi escrito em linguagem C com biblioteca de paralelização MPI.