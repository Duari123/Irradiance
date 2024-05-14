from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt

#Configurando fonte e tamanho da fonte dos plots
plt.rcParams.update({'font.family':'serif','font.size': 23})

plt.close('all') #Fecha figuras

# =============================================================================
# =============================================================================

def dia_do_ano():
    data_atual = datetime.now()
    dia_do_ano = (data_atual - datetime(data_atual.year, 1, 1)).days + 1
    """ Caso deseje utilizar o dia do ano correspondente a data atual, basta descomentar dia_do_ano no return """
    return 100#dia_do_ano

def hora_atual_fracionada():
    data_atual = datetime.now()
    hora_fracionada = data_atual.hour + data_atual.minute / 60 + data_atual.second / 3600
    return hora_fracionada

def declinacaoSolar(nda):
    delta = 23.45*np.sin( (360/365)*(nda - 81)*(np.pi/180) )
    return delta

def angulo_horario(hora_fracionada):
    angulo_h = 15 * (hora_fracionada - 12)
    return angulo_h

def converter_para_radianos(latitude, longitude):
    latitude_radianos = np.deg2rad(latitude)
    longitude_radianos = np.deg2rad(longitude)
    return latitude_radianos, longitude_radianos

def desvio_zenital(declinacao, latitude, angulo_horario):
    declinacao_radianos = np.radians(declinacao)
    latitude_radianos = np.radians(latitude)
    angulo_horario_radianos = np.radians(angulo_horario)

    cos_Z = np.sin(declinacao_radianos) * np.sin(latitude_radianos) + \
            np.cos(declinacao_radianos) * np.cos(latitude_radianos) * np.cos(angulo_horario_radianos)

    desvio_zenital_radianos = np.arccos(cos_Z)

    return desvio_zenital_radianos

def equacao_do_tempo(dia_do_ano):
    # B é uma aproximação da posição angular da Terra na sua órbita
    B = 360/365 * (dia_do_ano - 81)
    B_radianos = np.deg2rad(B)

    # Calculando a Equação do Tempo em minutos
    E = 9.87 * np.sin(2 * B_radianos) - 7.53 * np.cos(B_radianos) - 1.5 * np.sin(B_radianos)

    return E

def meio_dia_solar(dia_do_ano, longitude, fuso_horario):
    # Calculando a Equação do Tempo (E) em minutos
    E = equacao_do_tempo(dia_do_ano)

    # Correção de longitude em minutos
    correcao_longitude = 4 * (longitude - fuso_horario * 15)

    # Calculando o meio-dia solar local em minutos
    meio_dia_local_minutos = 12 * 60 - (E + correcao_longitude)

    # Convertendo minutos para horas fracionadas
    meio_dia_local_fracionado = meio_dia_local_minutos / 60

    return meio_dia_local_fracionado

def calcular_nascer_do_sol(meio_dia_solar, hsr):
    hsr_em_horas = hsr / 15
    nascer_do_sol_fracionado = meio_dia_solar - hsr_em_horas
    return nascer_do_sol_fracionado

def calcular_por_do_sol(meio_dia_solar, hsr):
    hsr_em_horas = hsr / 15
    por_do_sol_fracionado = meio_dia_solar + hsr_em_horas
    return por_do_sol_fracionado

def calculaLT(E, longitude):
    Lt = 60*( -E/60 - (4/60)*(45 - abs(longitude)) )
    return Lt

def calcular_hsr(latitude, declinacao):
    latitude_rad = np.radians(latitude)
    declinacao_rad = np.radians(declinacao)
    
    hsr_rad = np.arccos(-np.tan(latitude_rad) * np.tan(declinacao_rad))
    hsr_deg = np.degrees(hsr_rad)
    
    return hsr_deg

def hora_angular(hora_atual_fracionada, meio_dia_solar):
    return 15 * (hora_atual_fracionada - meio_dia_solar)

def calcular_beta(latitude, declinacao, hora_angular):
    latitude_rad = np.radians(latitude)
    declinacao_rad = np.radians(declinacao)
    hora_angular_rad = np.radians(hora_angular)

    beta_rad = np.arcsin(
        np.sin(declinacao_rad) * np.sin(latitude_rad)
        + np.cos(declinacao_rad) * np.cos(latitude_rad) * np.cos(hora_angular_rad)
    )

    return np.degrees(beta_rad)

def calcular_azimutal(declinacao, hora_angular, beta):
    declinacao_rad = np.radians(declinacao)
    hora_angular_rad = np.radians(hora_angular)
    beta_rad = np.radians(beta)

    sin_A = (np.cos(declinacao_rad) * np.sin(hora_angular_rad)) / np.cos(beta_rad)
    A_rad = np.arcsin(sin_A)

    # Convertendo de radianos para graus
    A = np.degrees(A_rad)

    # Ajustando o ângulo dependendo do quadrante
    if np.isscalar(hora_angular):
        if hora_angular > 0:
            A = 360 - A
    else:
        A[hora_angular > 0] = 360 - A[hora_angular > 0]

    return A


def coeficiente_backtracking(beta, distancia, largura):
    b = ( distancia/largura )*np.sin( np.deg2rad(beta) )
    return b

def calcula_gama(b):
    gama = 0
    if b <= 1 and b >= -1:
        gama = np.arcsin(b)
    elif b > 1:
        gama = 0
     
    return gama

def calcula_incidencia_solar(hora_atual_frac, hora_nascer_sol_frac, hora_por_sol_frac, HSR_graus):
    VD2284 = hora_atual_frac - hora_nascer_sol_frac
    VD2292 = 180*VD2284
    VD2296 = HSR_graus*(2/15)
    VD2300 = VD2292/VD2296
    IncidenciaSolar = VD2300 - 90
    
    if hora_atual_frac <= hora_nascer_sol_frac or hora_atual_frac >= hora_por_sol_frac:
        IncidenciaSolar = 0
    
    return IncidenciaSolar

def simula_pos(cof_back, hora_atual_frac, NascerSol, PorSol, meio_dia_solar, gama, incidencia_solar):
    ang_corrigido = 0
    ang_correcao = 0
    
    if cof_back <= 1 and hora_atual_frac <=  meio_dia_solar:
        ang_correcao = 90 - gama
    if cof_back <= 1 and hora_atual_frac >=  meio_dia_solar:
        ang_correcao = gama - 90
    if cof_back > 1:
        ang_correcao = 0
    
    ang_corrigido = ang_correcao + incidencia_solar
    
    if ang_corrigido >= 45:
        ang_corrigido = 45
    elif ang_corrigido <= -45:
        ang_corrigido = -45
    
    if hora_atual_frac <= NascerSol:
        ang_corrigido = 0
    elif hora_atual_frac >= PorSol:
        ang_corrigido = 0
    
    return ang_corrigido

# latitude = -19.9
# latitudeNp = np.array([-3.11, -10.11, -20.11, -30.15, -40])
latitudeNp = np.array([-20.14])
longitude = -44.8901
fuso_horario = -3 # Fuso horário de Brasília

for f in range(len(latitudeNp)):
    print("Número do dia do ano:", dia_do_ano())
    print("Hora atual fracionada:", hora_atual_fracionada())
    print("Declinação solar (graus):", declinacaoSolar(dia_do_ano()))
    print("Ângulo Horário Solar (graus):", angulo_horario(hora_atual_fracionada()))
    print("Desvio Zenital (graus):", np.rad2deg(desvio_zenital( declinacaoSolar(dia_do_ano()),  latitudeNp[f], angulo_horario(hora_atual_fracionada()))))
    oz = 1*(desvio_zenital( declinacaoSolar(dia_do_ano()),  latitudeNp[f], angulo_horario(hora_atual_fracionada())))
    print("Valor de 'E' (minutos):", equacao_do_tempo (dia_do_ano() ))
    print("Valor de 'LT':", calculaLT(equacao_do_tempo (dia_do_ano() ),longitude))
    print("Valor de 'HSR' (graus):", calcular_hsr(latitudeNp[f],declinacaoSolar(dia_do_ano())))
    
    meio_dia_fracionado = meio_dia_solar(dia_do_ano(), longitude, fuso_horario)
    print(f"Meio-dia solar local: {meio_dia_fracionado:.2f} horas")
    
    print("")
    
    print("Declinação solar (rad):", np.deg2rad( declinacaoSolar(dia_do_ano() ) ))
    print("Ângulo Horário Solar (rad):", np.deg2rad( angulo_horario(hora_atual_fracionada() ) ))
    print("Desvio Zenital (rad):", desvio_zenital( declinacaoSolar(dia_do_ano()),  latitudeNp[f], angulo_horario(hora_atual_fracionada())))
    
    print("")
    
    lat_radianos, lon_radianos = converter_para_radianos(latitudeNp[f], longitude)
    print(f"Latitude em radianos: {lat_radianos}")
    print(f"Longitude em radianos: {lon_radianos}")
    
    print("")
    
    # Exemplo de uso:
    meio_dia_fracionado = meio_dia_solar(dia_do_ano(), longitude, fuso_horario)
    hsr = calcular_hsr(latitudeNp[f], declinacaoSolar(dia_do_ano()))
    nascer_do_sol_fracionado = calcular_nascer_do_sol(meio_dia_fracionado, hsr)
    
    horas_nascer_do_sol = int(nascer_do_sol_fracionado)
    minutos_nascer_do_sol = int((nascer_do_sol_fracionado - horas_nascer_do_sol) * 60)
    
    print(f"Nascer do Sol (fracionado): {nascer_do_sol_fracionado:.2f} horas")
    print(f"Nascer do Sol (não fracionado): {horas_nascer_do_sol}:{minutos_nascer_do_sol}")
    
    print("")
    
    # Exemplo de uso:
    meio_dia_fracionado = meio_dia_solar(dia_do_ano(), longitude, fuso_horario)
    por_do_sol_fracionado = calcular_por_do_sol(meio_dia_fracionado, hsr)
    
    horas_por_do_sol = int(por_do_sol_fracionado)
    minutos_por_do_sol = int((por_do_sol_fracionado - horas_por_do_sol) * 60)
    
    print(f"Pôr do Sol (fracionado): {por_do_sol_fracionado:.2f} horas")
    print(f"Pôr do Sol (não fracionado): {horas_por_do_sol}:{minutos_por_do_sol}")
    
    print("")
    
    hora_atual = hora_atual_fracionada()
    ha = hora_angular(hora_atual, meio_dia_fracionado)
    ha_r = np.deg2rad(ha)
    
    print(f"Hora Angular: {ha:.2f} graus")
    print(f"Hora Angular: {ha_r:.2f} rad")
    
    print("")
    
    declinacao = declinacaoSolar(dia_do_ano())
    beta = calcular_beta(latitudeNp[f], declinacao, ha)
    beta_r = np.deg2rad(beta)
    
    print(f"Ângulo de elevação solar (Beta): {beta:.2f} graus")
    print(f"Ângulo de elevação solar (Beta): {beta_r:.2f} rad")
    
    print("")
    
    azimutal = calcular_azimutal(declinacao, ha, beta)
    
    print(f"Ângulo Azimutal: {azimutal:.2f} graus")
    
    print("")
    
    distancia = 6
    largura = 2.3
    b = coeficiente_backtracking(beta, distancia, largura)
    gama_r = calcula_gama(b) #Radianos
    gama = np.rad2deg(gama_r) #Graus
    
    print(f"Coeficiente de Backtracking: {b}")
    print(f"Gama (graus): {gama}")
    print(f"Gama (rad): {gama_r}")
    
    print("")
    
    HSR_graus = calcular_hsr(latitudeNp[f],declinacaoSolar(dia_do_ano()))
    IncSolar = calcula_incidencia_solar(hora_atual_fracionada(), nascer_do_sol_fracionado, por_do_sol_fracionado, HSR_graus)
    print(f"Incidência Solar: {IncSolar}")
    print("")
    
    posSimul = simula_pos(b, hora_atual_fracionada(), nascer_do_sol_fracionado, por_do_sol_fracionado ,meio_dia_fracionado, gama, IncSolar)
    print(f"Posição Simulada: {posSimul}")

# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================

# Função para calcular a posição do tracker em um momento específico
def calcular_posicao(hora_fracionada, latitude, longitude, fuso_horario, distancia, largura):
    declinacao = declinacaoSolar(dia_do_ano())
    hsr = calcular_hsr(latitude, declinacao)
    meio_dia_fracionado = meio_dia_solar(dia_do_ano(), longitude, fuso_horario)
    nascer_do_sol_fracionado = calcular_nascer_do_sol(meio_dia_fracionado, hsr)
    por_do_sol_fracionado = calcular_por_do_sol(meio_dia_fracionado, hsr)
    ha = hora_angular(hora_fracionada, meio_dia_fracionado)
    beta = calcular_beta(latitude, declinacao, ha)
    b = coeficiente_backtracking(beta, distancia, largura)
    gama = np.rad2deg(calcula_gama(b))
    IncSolar = calcula_incidencia_solar(hora_fracionada, nascer_do_sol_fracionado, por_do_sol_fracionado, hsr)
    posSimul = simula_pos(b, hora_fracionada, nascer_do_sol_fracionado, por_do_sol_fracionado, meio_dia_fracionado, gama, IncSolar)
    
    return posSimul

# Lista para armazenar os resultados
posicoes = []
Incidncia = []

# Lista de diferentes valores de tam que você deseja simular
# tams = [35, 40, 50, 60, 70, 80, 90, 100]  # Adicione ou modifique conforme necessário
tams = [50]
# Dicionários para armazenar os resultados
dict_posicoes = {}
dict_incidencias = {}

for lat in latitudeNp:
    # Recalcular as variáveis dependentes da latitude
    declinacao = declinacaoSolar(dia_do_ano())
    hsr = calcular_hsr(lat, declinacao)
    meio_dia_fracionado = meio_dia_solar(dia_do_ano(), longitude, fuso_horario)
    nascer_do_sol_fracionado = calcular_nascer_do_sol(meio_dia_fracionado, hsr)
    por_do_sol_fracionado = calcular_por_do_sol(meio_dia_fracionado, hsr)

    # Loop secundário para iterar sobre cada valor de tam
    for tam in tams:
        posicoes_temp = []
        incidencias_temp = []

        # Loop interno para simular a posição do rastreador solar a cada hora do dia
        for hora_fracionada in np.linspace(nascer_do_sol_fracionado, por_do_sol_fracionado, tam):
            posicao = calcular_posicao(hora_fracionada, lat, longitude, fuso_horario, distancia, largura)
            posicoes_temp.append(posicao)

            Incd = calcula_incidencia_solar(hora_fracionada, nascer_do_sol_fracionado, por_do_sol_fracionado, hsr)
            if Incd >= 45:
                Incd = 45
            if Incd <= -45:
                Incd = -45
            incidencias_temp.append(Incd)

        # Adicionando os resultados ao dicionário
        dict_posicoes[(lat, tam)] = posicoes_temp
        dict_incidencias[(lat, tam)] = incidencias_temp

        # Plotando o gráfico para a latitude atual e valor atual de tam
        plt.figure()
        plt.plot(np.linspace(nascer_do_sol_fracionado, por_do_sol_fracionado, tam), posicoes_temp, linewidth='4.5', color='red', label='Com Backtracking')
        plt.plot(np.linspace(nascer_do_sol_fracionado, por_do_sol_fracionado, tam), incidencias_temp, '--k', linewidth='4.5', label='Sem Backtracking')
        plt.xlabel('Hora do Dia [fracionada]')
        plt.ylabel('Posição do Rastreador [graus]')
        plt.title(f'Ângulo Seguidor Solar para Latitude {lat} e Tamanho {tam} do Vetor')
        plt.legend(loc='best', prop={'size': 18.0})
        plt.grid(True)

def simular_posicoes_e_incidencias(tams, latitude, dia_do_ano, longitude, fuso_horario, distancia, largura):
    declinacao = declinacaoSolar(dia_do_ano)
    hsr = calcular_hsr(latitude, declinacao)
    meio_dia_fracionado = meio_dia_solar(dia_do_ano, longitude, fuso_horario)
    nascer_do_sol_fracionado = calcular_nascer_do_sol(meio_dia_fracionado, hsr)
    por_do_sol_fracionado = calcular_por_do_sol(meio_dia_fracionado, hsr)

    resultados = {}

    for tam in tams:
        horas_fracionadas = np.linspace(nascer_do_sol_fracionado, por_do_sol_fracionado, tam)
        posicoes_temp = []
        incidencias_temp = []

        for hora_fracionada in horas_fracionadas:
            posicao = calcular_posicao(hora_fracionada, latitude, longitude, fuso_horario, distancia, largura)
            posicoes_temp.append(posicao)

            Incd = calcula_incidencia_solar(hora_fracionada, nascer_do_sol_fracionado, por_do_sol_fracionado, hsr)
            if Incd >= 45:
                Incd = 45
            elif Incd <= -45:
                Incd = -45
            incidencias_temp.append(Incd)

        resultados[tam] = {
            'horas_fracionadas': horas_fracionadas,
            'posicoes': posicoes_temp,
            'incidencias': incidencias_temp
        }

    return resultados


# # Para acessar os dados de uma latitude e tam específico:
# lat_terceira_posicao = latitudeNp[0]
# tam_especifico = 90  # Por exemplo
# posicoes_terceira_latitude = dict_posicoes[(lat_terceira_posicao, tam_especifico)]
# incidencias_terceira_latitude = dict_incidencias[(lat_terceira_posicao, tam_especifico)]