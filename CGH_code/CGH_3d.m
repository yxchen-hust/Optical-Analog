clc;
clearvars;
close all;

length = 512;
distance = 0.3;
N = 3;
delta = 0.15;
lambda = 638e-9;  % 波长(500nm)
dx = 8e-6; % 采样间隔8um
x = 0:length - 1;
X = repmat(x, length, 1) - length/2;
y = x';
Y = repmat(y, 1, length) - length/2;
X = X*dx;
Y = Y*dx;

picture = imread('dog1.jpg');
gray_picture = rgb2gray(picture);            % 读入初始图,作为初始迭代分布振幅 
sz = size(picture);

max_itera = 6000;
RMS_GS = zeros(max_itera, 1);                       % 计算GS算法均方根误差

max_RMS = 0.01;

image_hologram = zeros(length, length);
image_goal_zi = cell(N);
rgb = cell(N);
gray = cell(N);
Amplitude = cell(N);
z = zeros(N);

for i = 1:N
    rgb{i} = imread(['dog', num2str(i), '.jpg']);               % 读入初始图 作为初始迭代分布振幅 
    gray{i} = double(rgb2gray(rgb{i}));
    Amplitude{i} = imresize(gray{i}, [length,length]);             % 适应调制器分辨率
    Amplitude{i} = Amplitude{i}./(max(max(Amplitude{i})));     % 归一化
    z(i) = i * delta + distance;
    phase = 2*pi*rand(length,length);                  % 产生随机相位
    image_goal_zi{i} = Amplitude{i}.*exp(1i*phase);            % GS算法初始复振幅分布
end

for n = 1:max_itera    % 设置最大迭代次数
    fprintf("iteration times: %d\n", n);
    for i = 1:N
        % GS算法 正向传播计算
        % 计算菲涅尔衍射
        k = 2*pi / lambda;  % 波数
        z_phase = exp(1i * k * z(i)) / (1i * lambda * z(i));  % 衍射相位因子

        H = z_phase * exp(1i * k / (2 * z(i)) * (X.^2 + Y.^2));  % 菲涅尔正向传播核
        H_fft = fftshift(fft2(ifftshift(H)));  % 菲涅尔正向传播核的Fourier变换

        % 通过卷积计算正向衍射光场
        U_fft_zi = fftshift(fft2(ifftshift(image_goal_zi{i}))) .* H_fft;  % 频域乘积
        image_hologram_zi = fftshift(ifft2(ifftshift(U_fft_zi)));  % 返回空间域
        
        image_hologram = image_hologram + image_hologram_zi;  % 将各个物面衍射后得到的全息面叠加
    end

    for i = 1:N
        % GS算法 逆向传播计算
        k = 2*pi / lambda;  % 波数
        z_phase = exp(1i * k * z(i)) / (1i * lambda * z(i));  % 衍射相位因子

        H_inv = (exp(-1i * k * z(i)) / (1i * lambda * z(i))) * exp(-1i * k / (2 * z(i)) * (X.^2 + Y.^2));  % 菲涅尔逆向传播核
        H_inv_fft = fftshift(fft2(ifftshift(H_inv)));  % 菲涅尔逆向传播核的Fourier变换

        hologram_angle = angle(image_hologram);

        % 通过卷积计算逆向衍射光场
        U_inv_fft_zi = fftshift(fft2(ifftshift(exp(1i * hologram_angle)))) .* H_inv_fft;
        image_goal_zi{i} = fftshift(ifft2(ifftshift(U_inv_fft_zi)));      % 逆傅里叶变换返回空域

        g_er = abs(Amplitude{i}) - abs(image_goal_zi{i})/max(max(abs(image_goal_zi{i})));     % 计算误差矩阵
        RMS_GS(n) = RMS_GS(n) + sqrt(mean2((g_er.^2)));        % 计算均方根误差

        image_goal_zi{i} = abs(Amplitude{i}).*(image_goal_zi{i}./abs(image_goal_zi{i})); % 直接用初始振幅约束，不做改变
    end

    if RMS_GS(n) < max_RMS
            fprintf("total iteration times: %d\n", n);
            break;
    end
end

angle_hologram = angle(image_hologram);
angle_hologram = mod(angle_hologram,2*pi);
gray_angle = uint8((angle_hologram / (2*pi)) * 255);
imwrite(gray_angle, '3d_cgh.png');

figure(1)
% 计算最终全息图生成的像 并与原图进行对比
for j = 1:N
    subplot(2, 3, j);
    imshow(Amplitude{j}, []);
    title(sprintf('原图%d（dog%d.jpg）', j, j));

    subplot(2, 3, j+N);
    imshow(abs(image_goal_zi{j}), []);
    title(sprintf('模拟衍射输出-第%d张图', j));
end

figure(2)
imshow(gray_angle,[]);title('相位原件分布(GS)');

figure(3)
plot(1:max_itera ,RMS_GS);
xlabel('循环次数');
ylabel('RMS误差(GS)');
