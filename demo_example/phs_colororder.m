function ColorOrder = pcag_colororder
% get colororder for figures

% ideas for colors: http://www.rapidtables.com/web/color/RGB_Color.htm

% RGB = [
%     0 0 204;
%     0 102 204;
%     0 204 204;
%     96 96 96;
%     204 204 0;
%     204 102 0;
%     204 0 0];

RGB = [
    0 0 153;
    0 0 255;
    51 153 255;
    0 0 0; % bar plot all: 40 40 41; bar plot 50%: 247 175 175; bar plot line: 148 151 153
    255 102 102;
    255 0 0;
    204 0 0];

ColorOrder = RGB./255; % convert to 0-1 instead of 0-256