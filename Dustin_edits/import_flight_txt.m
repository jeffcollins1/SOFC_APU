function [segment,history,profile] = import_flight_txt(filename)
if isempty(strfind(filename,'.txt'))
    filename = strcat(filename,'.txt');
end
text = fileread(filename);
n_line = strfind(text,char(10));
n_l = length(n_line);
lines = cell(n_l,1);
s = 1;
for i = 1:1:length(n_line)
    lines(i) = {text(s:n_line(i)-2)};
    s = n_line(i)+1;
end

i = 1;
j = 0;
while isempty(strfind(lines{i},'Manoeuvre segment'))
    i = i+1;
end
prop = {'initial_mass';'initial_alt';'delta_isa';'airspeed';'thrust_per_eng';'climb_rate';'flight_angle';'true_airspeed';'fuel_flow';'NOx';'HC';'CO';'LD_ratio';'total_drag';};
index = [1;2;3;4;8;10;11;12;14;15;16;17;19;20;];
hist = {'alt';'time';'dist';'burnt';'FN_eng';'cas';'mach';'RoC';};
while i<n_l-1
    if ~isempty(strfind(lines{i},'Manoeuvre segment'))
        j = j+1;
        i = i+1;
    end
    for k = 1:1:length(index)
        per = strfind(lines{i+index(k)},'.');
        val = lines{i+index(k)}(per(1)-7:per(1)+1);
        val = {strrep(val,'e','')};
        segment.(prop{k})(j,1) = str2double(val);
    end
    per = strfind(lines{i+5},'.');
    segment.flaps(j,1) = str2double(lines{i+5}(per(1)-9:per(1)-4));
    i = i+22;
    a = strfind(lines{i},'time = '); b = strfind(lines{i},'mins,');
    segment.end_time(j,1) = str2double(lines{i}(a(1)+7:b(1)-1));
    a = strfind(lines{i},'endmass = '); b = strfind(lines{i},'kg.');
    segment.end_mass(j,1) = str2double(lines{i}(a(1)+9:b(1)-1));
    i = i+9;
    t = 0;
    while isempty(strfind(lines{i},'---'))
        t = t+1;
        per = strfind(lines{i},'.');
        a = 1;
        for k = 1:1:length(hist)
            a = max(a,per(k)-6);
            b = min(length(lines{i}),per(k)+2);
            history.(hist{k})(j,t) = str2double(lines{i}(a:b));
            a = per(k)+1;
        end
        i = i+1;
    end
    while i<n_l-1 && isempty(strfind(lines{i},'Manoeuvre segment'))
        i = i+1;
    end
end
for k = 1:1:length(hist)
    profile.(hist{k}) = [];
end
for j = 1:1:length(history.alt(:,1))
    n = nnz(history.alt(j,:));
    for k = 1:1:length(hist)
        profile.(hist{k})(end+1:end+n,1) = history.(hist{k})(j,1:n)';
    end
end
figure(1)
[ax,h1,h2] = plotyy(profile.time/3600,max(0,profile.FN_eng),profile.time/3600,profile.alt);
xlabel(ax(1),'Time (hrs)')
ylabel(ax(1),'Thrust per Engine (N)')
ylabel(ax(2),'Altitude (ft)')
end%Ends function import_flight_txt