# aerodynamic_moment_model2_2_MCoeff_optimal

## 注意事项：该文档针对翅膀扭转轴向前缘移动之后

 aerodynamic_moment_model2_2_optimal——用于气动力系数的优化——四个变量

 (1)正确的无量纲弦长计算——计算正确;

 (2)正确的有效力臂的无量纲位置Z_nd计算——计算正确;

 (3)在扭转轴向上偏移C_maxy之后——用最大前缘点和最小后缘点弦向坐标差对应的片条长度来计算扭转轴的位置 ——在最大前缘后面的0.25*[最大前缘点y坐标和最小后缘点y坐标(即弦向坐标差)=C_max]位置时;

 (4) 翅平面坐标系的原点与扭转轴的近端点重合.

## 1.计算翅膀形貌学参数;

## 2.计算气动力和力矩，以及刚性翅转动时的惯性力和离心力.

![aerodynamic_moment](https://github.com/xijunke/aerodynamic_moment_model2_2_MCoeff_optimal/blob/master/pic_png_tif_eps_pdf/%E8%AE%A1%E7%AE%97%E5%92%8C%E5%AE%9E%E6%B5%8B%E8%99%AB%E4%BD%93%E5%9D%90%E6%A0%87%E4%B8%8B%E7%9A%84%E4%BF%AF%E4%BB%B0%E5%8A%9B%E7%9F%A9%E7%9A%84%E5%AF%B9%E6%AF%943.png)

<div align=center>
<img src="https://github.com/xijunke/aerodynamic_moment_model2_2_MCoeff_optimal/blob/master/pic_png_tif_eps_pdf/%E6%B0%94%E5%8A%A8%E5%8A%9B%E7%9F%A9%E5%88%86%E9%87%8F%E9%9A%8F%E6%97%B6%E9%97%B4%E7%9A%84%E5%8F%98%E5%8C%96%E8%A7%84%E5%BE%8B%E4%B8%8E%E5%AE%9E%E9%AA%8C%E6%B5%8B%E8%AF%95%E7%BB%93%E6%9E%9C%E7%9A%84%E5%AF%B9%E6%AF%94%E2%80%94%E2%80%94%E4%BC%98%E5%8C%96%E4%B9%8B%E5%90%8E2.png" width="1200" height="700"/>
</div>