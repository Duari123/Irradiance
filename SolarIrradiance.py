# =============================================================================
# Importando Bibliotecas e código do Backtracking
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
import CodBackTracking as cbt

# =============================================================================
# Fechando todas as janelas de plot abertas
# =============================================================================

plt.close('all') #comando para fechar todas janelas de plot

# =============================================================================
# Configuração de layout
# =============================================================================

# Configurando fonte e tamanho da fonte dos plots
plt.rcParams.update({'font.family': 'serif', 'font.size': 21})

# =============================================================================
# Curva da Incidencia Solar
# =============================================================================

# Parâmetros de exemplo
tams = [50]  # Pode ser uma lista de diferentes tamanhos para simular
latitude = -22.14  # Exemplo de latitude
dia_do_ano = cbt.dia_do_ano()  # Exemplo de dia do ano
longitude = -44.8901  # Exemplo de longitude
fuso_horario = -3  # Exemplo de fuso horário
distancia = 6  # Distância entre os painéis
largura = 2.3  # Largura do painel

# Chama a função com os parâmetros desejados
resultados = cbt.simular_posicoes_e_incidencias(tams, latitude, dia_do_ano, longitude, fuso_horario, distancia, largura)

# Itera sobre os resultados para cada valor de tam
for tam, dados in resultados.items():
    horas_fracionadas = dados['horas_fracionadas']
    posicoes_temp = dados['posicoes']
    incidencias_temp = dados['incidencias']
    
# =============================================================================
# =============================================================================

def irradFixo(DiaAno, Latitude):
    L = Latitude
    n = DiaAno
    # Declinação solar
    delta = 23.45 * np.sin(np.radians((n - 81) * 360 / 365))

    # Ângulo horário
    Hsr = np.degrees(np.arccos(-1 * np.tan(np.radians(L)) * np.tan(np.radians(delta))))
    h_sunrise = 12 - Hsr / 15
    h_sunset = 12 + Hsr / 15
    h = np.arange(h_sunrise, h_sunset, 1 / (5 * 60))  # Passo de amostragem
    H = 15 * (h - 12)

    # Ângulo de altitude do sol
    beta = np.degrees(np.arcsin(np.cos(np.radians(L)) * np.cos(np.radians(delta)) * np.cos(np.radians(H)) + np.sin(np.radians(L)) * np.sin(np.radians(delta))))
    seno_azimuth = (np.cos(np.radians(delta)) * np.sin(np.radians(H))) / np.cos(np.radians(beta))
    fi_azimuth = np.degrees(np.arcsin(seno_azimuth))

    # Lógica para ajustar azimuth
    tt = np.append([0.2], np.diff(fi_azimuth)) / (1 / (5 * 60))
    crit = round(len(fi_azimuth)) / 2
    cont = 0
    if max(seno_azimuth) >= 0.9999:
        for j in range(len(tt)):
            cont += 1
            if tt[j] <= 0 and cont < crit:
                fi_azimuth[j] = -1 * fi_azimuth[j] - 180
            if tt[j] <= 0 and cont > crit:
                fi_azimuth[j] = -1 * fi_azimuth[j] + 180

    # Curvas de irradiância
    A = 1160 + 75 * np.sin(np.radians((n - 275) * 360 / 365))
    k = 0.174 + 0.035 * np.sin(np.radians((n - 100) * 360 / 365))
    m = 1 / np.sin(np.radians(beta))
    Ib = A * np.exp(-k * m)
    C = 0.095 + 0.04 * np.sin(np.radians((n - 100) * 360 / 365))
    Ib[np.isinf(Ib)] = 0

    # Painel Fixo
    Sigma = abs(round(L))
    theta = np.degrees(np.arccos(np.cos(np.radians(beta)) * np.cos(np.radians(fi_azimuth)) * np.sin(np.radians(Sigma)) + np.sin(np.radians(beta)) * np.cos(np.radians(Sigma))))
    IBC = Ib * np.cos(np.radians(theta))
    IDC = C * Ib * (1 + np.cos(np.radians(Sigma))) / 2
    IRC = ((1 - np.cos(np.radians(Sigma))) / 2) * 0.2 * (np.sin(np.radians(beta)) + C) * Ib
    Irr_fixo = IBC + IDC + IRC
    
    Irr_fixo = np.where(Irr_fixo < 0, 0, Irr_fixo)

    
    return Irr_fixo, h

# =============================================================================
# =============================================================================

def irrad1eixoCB_SB(DiaAno, LatitudePos, Tam, Modo):
    
    n = DiaAno
    #Passo da simulação -> Mantive o mesmo do Gusman
    step = 1 / (5 * 60)
    
    #Latitude local -> Alterado para Divinópolis - MG
    L = LatitudePos
    Tamanho = Tam
    
    #Reflexão devido ao solo
    rho = 0.2
    
    #Declinação solar -> Troquei para a função que criei previamente
    delta = cbt.declinacaoSolar(DiaAno)
    
    #Hora Angular -> Ao invés de deixar o meio-dia fixo em 12h, coloquei minha
    #função que o calcula dinamicamente
    Hsr = np.degrees(np.arccos(-1 * np.tan(np.radians(L)) * np.tan(np.radians(delta))))
    h_sunrise = cbt.meio_dia_fracionado - Hsr / 15
    h_sunset = cbt.meio_dia_fracionado + Hsr / 15
    h = np.arange(h_sunrise, h_sunset, step)
    H = 15 * (h - cbt.meio_dia_fracionado)
    
    #Variável Beta -> Troquei para a função que criei previamente    
    beta = cbt.calcular_beta(L, delta, H)
    
    #Achando Azimute
    seno_azimuth = (np.cos(np.radians(delta)) * np.sin(np.radians(H))) / np.cos(np.radians(beta))
    fi_azimuth = np.degrees(np.arcsin(seno_azimuth))
    
    #Ajustando quadrante do ângulo de Azimute
    tt = np.concatenate(([0.2], np.diff(fi_azimuth))) / step
    cont = 0
    crit = len(fi_azimuth) // 2
    if max(seno_azimuth) >= 0.9999:
        for j in range(len(tt)):
            cont += 1
            if tt[j] <= 0 and cont < crit:
                fi_azimuth[j] = -1 * fi_azimuth[j] - 180
            elif tt[j] <= 0 and cont > crit:
                fi_azimuth[j] = -1 * fi_azimuth[j] + 180
    
    #Irradiância extraterrestre com atenuação:**
    A = 1160 + 75 * np.sin((n - 275) * 2 * np.pi / 365)
    k = 0.174 + 0.035 * np.sin((n - 100) * 2 * np.pi / 365)
    m = 1 / np.sin(np.radians(beta))
    Ib = A * np.exp(-k * m)
    C = 0.095 + 0.04 * np.sin((n - 100) * 2 * np.pi / 365)
    Ib[np.isinf(Ib)] = 0
    
    if Modo == 0:
    
        Sigma = []
        X = len(H)
        
        # Chama a função com os parâmetros desejados
        resultados = cbt.simular_posicoes_e_incidencias(Tamanho, L, DiaAno, longitude, fuso_horario, distancia, largura)

        # Itera sobre os resultados para cada valor de tam
        for tam, dados in resultados.items():
            posicoes_temp = dados['posicoes']
            
        V = np.array( posicoes_temp )*-1
        
        for v in V:
            Sigma.extend([v] * (X // len(V)))
        while len(Sigma) < X:
            Sigma.append(45)
        theta = np.degrees(np.arccos(np.cos(np.radians(beta)) * np.cos(np.radians(fi_azimuth + 90)) * np.sin(np.radians(Sigma)) + np.sin(np.radians(beta)) * np.cos(np.radians(Sigma))))
        IBC_1x = Ib * np.cos(np.radians(theta))
        IDC_1x = C * Ib * ((1 + np.cos(np.radians(90 - beta + delta))) / 2)
        IRC_1x = ((1 - np.cos(np.radians(90 - beta + delta))) / 2) * rho * (np.sin(np.radians(beta)) + C) * Ib
        Irr_1eixoComBackTracking = IBC_1x + IDC_1x + IRC_1x
        
        for i in range(len(Irr_1eixoComBackTracking)):
            if Irr_1eixoComBackTracking[i] < 0:
                Irr_1eixoComBackTracking[i] = 0
    
        return Irr_1eixoComBackTracking, h
    
    elif Modo == 1:
        Sigma = []
        X = len(H)
        
        # Chama a função com os parâmetros desejados
        resultados = cbt.simular_posicoes_e_incidencias(Tamanho, L, DiaAno, longitude, fuso_horario, distancia, largura)

        # Itera sobre os resultados para cada valor de tam
        for tam, dados in resultados.items():
            posicoes_temp = dados['posicoes']
            incidencias_temp = dados['incidencias']
            
        S = 45 
        V = np.array( incidencias_temp )*-1
        
        for v in V:
            Sigma.extend([v] * (X // len(V)))
        while len(Sigma) < X:
            Sigma.append(-S)
        theta = np.degrees(np.arccos(np.cos(np.radians(beta)) * np.cos(np.radians(fi_azimuth + 90)) * np.sin(np.radians(Sigma)) + np.sin(np.radians(beta)) * np.cos(np.radians(Sigma))))
        IBC_1x = Ib * np.cos(np.radians(theta))
        IDC_1x = C * Ib * ((1 + np.cos(np.radians(90 - beta + delta))) / 2)
        IRC_1x = ((1 - np.cos(np.radians(90 - beta + delta))) / 2) * rho * (np.sin(np.radians(beta)) + C) * Ib
        Irr_1eixoSemBackTracking = IBC_1x + IDC_1x + IRC_1x
        
        for i in range(len(Irr_1eixoSemBackTracking)):
            if Irr_1eixoSemBackTracking[i] < 0:
                Irr_1eixoSemBackTracking[i] = 0
                
        Irr_1eixoSemBackTracking = np.where(Irr_1eixoSemBackTracking < 0, 0, Irr_1eixoSemBackTracking)
        
        return Irr_1eixoSemBackTracking, h

# =============================================================================
# =============================================================================

def irrad1eixo(DiaAno, Latitude, PassosTamTrack):
    n = DiaAno
    step = 1 / (5 * 60)  # Passo da simulação
    
    L = Latitude  # Latitude local

    # Declinação solar
    delta = 23.45 * np.sin(np.radians((n - 81) * 360 / 365))
    
    # Hora Angular
    Hsr = np.degrees(np.arccos(-np.tan(np.radians(L)) * np.tan(np.radians(delta))))
    h_sunrise = cbt.meio_dia_fracionado - Hsr / 15
    h_sunset = cbt.meio_dia_fracionado + Hsr / 15
    h = np.arange(h_sunrise, h_sunset, step)
    H = 15 * (h - cbt.meio_dia_fracionado)

    # Variável Beta
    beta = np.degrees(np.arcsin(np.cos(np.radians(L)) * np.cos(np.radians(delta)) * np.cos(np.radians(H)) + np.sin(np.radians(L)) * np.sin(np.radians(delta))))
    
    # Achando Azimute
    seno_azimuth = (np.cos(np.radians(delta)) * np.sin(np.radians(H))) / np.cos(np.radians(beta))
    fi_azimuth = np.degrees(np.arcsin(seno_azimuth))

    # Ajustando quadrante do ângulo de Azimute
    tt = np.concatenate(([0.2], np.diff(fi_azimuth))) / step
    for j in range(len(tt)):
        if tt[j] <= 0:
            fi_azimuth[j] = 180 - fi_azimuth[j] if j < len(fi_azimuth) // 2 else -180 - fi_azimuth[j]
    
    # Irradiância extraterrestre com atenuação
    A = 1160 + 75 * np.sin(np.radians((n - 275) * 360 / 365))
    k = 0.174 + 0.035 * np.sin(np.radians((n - 100) * 360 / 365))
    C = 0.095 + 0.05 * np.sin(np.radians((n - 100) * 360 / 365))
    m = 1 / np.sin(np.radians(beta))
    Ib = A * np.exp(-k * m)
    Ib[np.isinf(Ib)] = 0

    passo_tracker = [PassosTamTrack]
    for jj in passo_tracker:
        X = len(H)
        S = 45  # Ângulo máximo
        V = np.linspace(S, -S, jj)
        Sigma = np.repeat(V, np.floor(X / jj))
        Sigma = np.pad(Sigma, (0, X - len(Sigma)), 'constant', constant_values=(-S))
    
        Sig = np.sort(H)[::-1]
        theta = np.degrees(np.arccos(np.cos(np.radians(beta)) * np.cos(np.radians(fi_azimuth + 90)) * np.sin(np.radians(Sig)) + np.sin(np.radians(beta)) * np.cos(np.radians(Sig))))
        IBC_1x = Ib * np.cos(np.radians(theta))
        IDC_1x = C*Ib*((1+np.cos(np.radians(90-beta+delta)))/2);
        IRC_1x = ((1-np.cos(np.radians(90-beta+delta)))/2)*0.2*(np.sin(np.radians(beta))+C)*Ib;
        Irr_1eixoMax = IBC_1x + IDC_1x + IRC_1x
    
    Irr_1eixoMax = np.where(Irr_1eixoMax < 0, 0, Irr_1eixoMax)
    
    passo_tracker = [PassosTamTrack]
    Irr_1eixoMPassos = None
    for jj in passo_tracker:
        X = len(H)
        S = 45
        V = np.linspace(S, -S, jj)
        Sigma = np.repeat(V, np.floor(X / jj))
        Sigma = np.pad(Sigma, (0, X - len(Sigma)), 'constant', constant_values=(-S))

        theta = np.degrees(np.arccos(np.cos(np.radians(beta)) * np.cos(np.radians(fi_azimuth + 90)) * np.sin(np.radians(Sigma)) + np.sin(np.radians(beta)) * np.cos(np.radians(Sigma))))
        IBC_1x = Ib * np.cos(np.radians(theta))
        IDC_1x = C * Ib * ((1 + np.cos(np.radians(90 - beta + delta))) / 2)
        IRC_1x = ((1 - np.cos(np.radians(90 - beta + delta))) / 2) * 0.2 * (np.sin(np.radians(beta)) + C) * Ib
        Irr_1eixo = IBC_1x + IDC_1x + IRC_1x

        if jj == max(passo_tracker):
            Irr_1eixoMPassos = Irr_1eixo

    return Irr_1eixoMax, h, Irr_1eixoMPassos

# =============================================================================
# Função para Calcular Irradiância Tracker 2 Eixos
# =============================================================================

def irrad2eixos(DiaAno, Latitude):
    
    n = DiaAno
    #Passo da simulação -> Mantive o mesmo do Gusman
    step = 1 / (5 * 60)
    
    #Latitude local -> Alterado para Divinópolis - MG
    L = Latitude
    
    #Reflexão devido ao solo
    rho = 0
    
    #Declinação solar -> Troquei para a função que criei previamente
    delta = cbt.declinacaoSolar(DiaAno)
    
    #Hora Angular -> Ao invés de deixar o meio-dia fixo em 12h, coloquei minha
    #função que o calcula dinamicamente
    Hsr = np.degrees(np.arccos(-1 * np.tan(np.radians(L)) * np.tan(np.radians(delta))))
    h_sunrise = cbt.meio_dia_fracionado - Hsr / 15
    h_sunset = cbt.meio_dia_fracionado + Hsr / 15
    h = np.arange(h_sunrise, h_sunset, step)
    H = 15 * (h - cbt.meio_dia_fracionado)
    
    #Variável Beta -> Troquei para a função que criei previamente    
    beta = cbt.calcular_beta(L, delta, H)
    
    #Achando Azimute
    seno_azimuth = (np.cos(np.radians(delta)) * np.sin(np.radians(H))) / np.cos(np.radians(beta))
    fi_azimuth = np.degrees(np.arcsin(seno_azimuth))
    
    # Ajustando quadrante do ângulo de Azimute
    tt = np.concatenate(([0.2], np.diff(fi_azimuth))) / step
    for j in range(len(tt)):
        if tt[j] <= 0:
            fi_azimuth[j] = 180 - fi_azimuth[j] if j < len(fi_azimuth) // 2 else -180 - fi_azimuth[j]
    
    # Irradiância extraterrestre com atenuação:
    A = 1160 + 75 * np.sin(np.radians((n - 275) * 360 / 365))
    k = 0.174 + 0.035 * np.sin(np.radians((n - 100) * 360 / 365))
    m = 1 / np.sin(np.radians(beta))
    Ib = A * np.exp(-k * m)
    Ib[np.isinf(Ib)] = 0
    C = 0.095 + 0.05 * np.sin(np.radians((n - 100) * 360 / 365))
    
    passo_tracker = [10]
    for jj in passo_tracker:
        X = len(H)
        S = 45  # Ângulo máximo
        V = np.linspace(S, -S, jj)
        Sigma = np.repeat(V, np.floor(X / jj))
        Sigma = np.pad(Sigma, (0, X - len(Sigma)), 'constant', constant_values=(-S))

        IBC_2x = Ib
        IDC_2x = C * Ib * ((1 + np.cos(np.radians(90 - beta))) / 2)
        IRC_2x = ((1 - np.cos(np.radians(90 - beta))) / 2) * rho * (np.sin(np.radians(beta)) + C) * Ib
        Irr_2eixoMax = IBC_2x + IDC_2x + IRC_2x
    
    Irr_2eixoMax = np.where(Irr_2eixoMax < 0, 0, Irr_2eixoMax)
    
    return Irr_2eixoMax, h

# =============================================================================
# Plot's de exemplo
# =============================================================================

"""
Caso fixo, rastreadores de um e dois eixos (irradiâncias máximas para um dia)
"""

# Calculando irradiância para a configuração fixa
Irr_fixo, h_fixo = irradFixo(dia_do_ano, latitude)
Irr_1eixo, h_1eixo, Irr_1eixoSteps = irrad1eixo(dia_do_ano, latitude, 10)
Irr_2eixo, h_2eixo = irrad2eixos(dia_do_ano, latitude)

# Plotagem dos resultados
plt.figure(figsize=(10, 6))
plt.plot(h_fixo, Irr_fixo, label='Sistema Fixo', color='blue')
plt.plot(h_1eixo, Irr_1eixo, label='Tracker 1-eixo', color='black')
plt.plot(h_2eixo, Irr_2eixo, label='Tracker 2-eixos', color='red')
plt.xlabel('Hora do dia')
plt.ylabel('Irradiância (W/m²)')
plt.title(f'Irradiância Solar ao longo do dia - Lat: {latitude}, Dia: {dia_do_ano}')
plt.legend()

"""
Rastreador de um eixo para diferentes categorias de movimentação
"""

# Passos Fixos, 10 ao longo do dia
Irr_1eixo, h_1eixo, Irr_1eixoSteps = irrad1eixo(dia_do_ano, latitude, 10)

# Plotagem dos resultados
plt.figure(figsize=(10, 6))
plt.plot(h_1eixo, Irr_1eixo, label='Irrad. Max.', color='black')
plt.plot(h_fixo, Irr_1eixoSteps, label='Irrad. Com Passos Fixos', color='blue')
plt.xlabel('Hora do dia')
plt.ylabel('Irradiância (W/m²)')
plt.title(f'Irradiância Solar ao longo do dia para tracker solar de um eixo - Lat: {latitude}, Dia: {dia_do_ano}')
plt.legend()

# Movimento Adptivo - Sem Backtracking
Steps = 55
Irr_1eixoSemBackTracking, hsb = irrad1eixoCB_SB(dia_do_ano, latitude, [Steps], 1)
Irr_1eixoComBackTracking, hcb = irrad1eixoCB_SB(dia_do_ano, latitude, [Steps], 0)

# Plotando o gráfico
plt.figure(figsize=(10, 6))
plt.title(f"Irradiância Tracker Mov. Adptativo [Lat.: {latitude}; Tam.: {Steps}; Dia: {dia_do_ano}]")
plt.plot(hsb, Irr_1eixoSemBackTracking, 'k', linewidth = '5.5', label='Sem Backtracking')
plt.plot(hcb, Irr_1eixoComBackTracking, 'r', linewidth = '5.5', label='Com Backtracking')
plt.ylabel('Irradiância (W/m²)')
plt.xlabel('Hora do dia')
plt.legend()