clear all;
T = 8; % Часы моделирования

nf = 2; % количество факторов в плане эксперимента
minf = [4 0.3]; % минимальное значение факторов a и b
maxf = [16 2.4]; % максимальное значение факторов a и b

% Стратегическое планирование эксперимента
% Формирование дробного двухуровневого плана эксперимента для учета взаимодействий
fracfact('a b ab'); % формируем список возможных взаимодействий
% вычисление общего количества точек в плане эксперимента 
% (все возможные сочетания уровней факторов)
N = 2 ^ nf;
fracplan = ans;

% Сформируем матрицу X с добавлением столбца значений фиктивного фактора
fictfact = ones(N, 1); % ones - возвращает скалярную 1, либо матрицу этих единиц
X = [fictfact ans]'; % формируем транспонированную матрицу плана с добавлением фиктивного фактора

fraceks = zeros(N, nf); % матрица (кол-во экспериментов x кол-во факторов)
for i = 1:nf
    for j = 1:N
        % Заполняем матрицу значениями
        fraceks(j,i) = minf(i) + (fracplan(j,i) + 1) * (maxf(i) - minf(i)) / 2;
    end
end

p0 = zeros(1, N); % относительная пропускная способность
A = zeros(1, N); % абсолютная пропускная способность
p_refuse = zeros(1, N); % вероятность отказа обслуживания

% Сколько телефонов должно быть в агентстве, чтобы относительная 
% пропускная способность была не менее 0,75
n = zeros(1, N);

for i= 1:N
    % Интенсивность в этот час
    lambda = fraceks(i, 1);
    t = fraceks(i, 2);
    % Моделирование
    s = sim('trenl', 60*T);
    % интенсивность входящего потока, заявок в минуту
    average_calls_per_minute = s.calls / 60;
    % интенсивность потока обслуживания, заявок за минуту
    average_servs_per_minute = 1 / s.average_serv_time;
    % Относительная пропускная способность (телефонная линия свободна, заявок нет)
    p0(i) = average_servs_per_minute / (average_servs_per_minute + average_calls_per_minute);
    % Абсолютная пропускная способность (заявок, обслуживаемых в минуту)
    A(i) = p0(i) * average_calls_per_minute;
    % вероятность отказа (занятости телефона)
    p_refuse(i) = 1 - p0(i);
    % Сколько телефонов должно быть в агентстве, чтобы относительная 
    % пропускная способность была не менее 0,75
   
    % приведенная интенсивность входящего потока
    ro = average_calls_per_minute / average_servs_per_minute;
    % cdf
    cdf = 0;
    % количество обслуживающих телефонов
    k = 0;
    while(0.75 > cdf)
        k = k + 1;
        p = 0;
        for j=0:k
            p = p + ((ro ^ j) / factorial(j));
        end
        p = 1 / p;
        % Вероятность отказа k-го телефона
        p_k = p * ((ro ^ k) / factorial(k));
        % пропускная способность системы с k-количестом телефонов
        cdf = 1 - p_k;
    end
    n(i) = k;
end

% Регрессионный анализ

% Определение коэффициентов регрессии
Cl = X * X';
b_l = inv(Cl) * X * p0';

% Формирование зависимости реакции системы на множестве реальных значений факторов
Al = minf(1):0.1:maxf(1);
Bl = minf(2):0.01:maxf(2);
[unused, N1] = size(Al);
[unused, N2] = size(Bl);

Yc = zeros(N2, N1);
for i = 1:N1
    for j = 1:N2
        anl = 2 * (Al(i) - minf(1)) / (maxf(1) - minf(1)) - 1;
        bnl = 2 * (Bl(j) - minf(2)) / (maxf(2) - minf(2)) - 1;
        % Экспериментальная поверхность реакции (линейная регрессия)
        Yc(j, i) = b_l(1) + anl * b_l(2) + bnl * b_l(3) + anl * bnl * b_l(4);
        % Теоритическая
        Yo(j, i) = (1 / Bl(j)) / ((Al(i) / 60) + (1 / Bl(j)));
    end
end

% Отображение зависимостей в трехмерной графике
[x, y] = meshgrid(Al, Bl); % координата в двумерном пространстве
figure;
% Исход в эксперименте
% Экспериментальная поверхность реакции
subplot(1, 2, 1), plot3(x, y, Yc);
xlabel('интенсивность N в час');
xlim([0 maxf(1) + 1]);
ylabel('время t');
ylim([0 maxf(2) + 0.1])
zlabel('Относительная пропускная способность');
zlim([0 1]);
title('Экспериментальная модель');
grid on;

% Теоретическая поверхность реакции
subplot(1, 2, 2), plot3(x, y, Yo);
xlabel('интенсивность N в час');
xlim([0 max(maxf(1)) + 1]);
ylabel('время t');
ylim([0 max(maxf(2)) + 0.1])
zlabel('Относительная пропускная способность');
zlim([0 1]);
title('Теоретическая модель');
grid on;

clear Al anl ans Bl bnl Cl fictfact fracplan i j lambda unused N N1 N2 nf t X ...
    k p p_k cdf ro current_lambda i j average_calls_per_minute ...
    average_servs_per_minute;