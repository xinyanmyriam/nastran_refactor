# NASTRAN-95 源码分类（基于 Zienkiewicz 有限元方法论）

## 分类依据

本分类基于 Zienkiewicz, Taylor & Zhu《The Finite Element Method: Its Basis and Fundamentals》
（第7版/第8版）的标准有限元计算流程，将 NASTRAN-95 `mis/` 文件夹中的子程序映射到
有限元方法的标准计算阶段。

**参考文献**：
- Zienkiewicz OC, Taylor RL, Zhu JZ. The Finite Element Method: Its Basis and Fundamentals. 7th ed. Butterworth-Heinemann, 2013.
- NASA SP-222. NASTRAN Programmer's Manual, 1970.

## Zienkiewicz 有限元计算标准流程

```
1. 问题离散化（Discretization）
   → 网格划分、节点编号
   
2. 单元公式化（Element Formulation）
   → 形函数、应变-位移关系、本构关系
   → 单元刚度矩阵 Ke = ∫ BᵀDB dΩ
   
3. 全局装配（Assembly）
   → 单元矩阵 → 全局矩阵 K, M, C
   
4. 边界条件施加（Boundary Conditions）
   → 约束处理、载荷向量
   
5. 方程组求解（Solution）
   → 静力学: Ku = F
   → 特征值: Kφ = λMφ
   → 瞬态: Mü + Cu̇ + Ku = F(t)
   
6. 后处理（Post-Processing）
   → 应力/应变恢复: σ = DBu
   → 结果输出
```

---

## NASTRAN-95 子程序分类

### 类别 A：单元刚度/质量矩阵生成（对应 Zienkiewicz Ch.5-9）

**A1. 1D 杆/桁架单元（Rod/Truss）**
| 文件 | 行数 | 功能 |
|------|:---:|------|
| krod.f | 308 | CROD 杆单元刚度矩阵（轴力+扭转）|
| ktube.f | 69 | CTUBE 管单元刚度矩阵 |

**A2. 1D 梁单元（Beam/Frame）**
| 文件 | 行数 | 功能 |
|------|:---:|------|
| kbar.f | 621 | CBAR 梁单元刚度矩阵（弯曲+轴力+扭转+剪切）|
| kelbow.f | 525 | CELBOW 弯管单元刚度矩阵 |

**A3. 2D 膜单元（Membrane/Plane Stress）**
| 文件 | 行数 | 功能 |
|------|:---:|------|
| ktrmem.f | 406 | CTRIA3 三角形膜单元刚度 |
| ktriqd.f | 210 | CQUAD4 四边形膜单元（驱动程序）|
| kqdmem.f | 288 | CQUAD4 四边形膜单元（核心计算）|
| ktrm6d.f | ? | 6节点三角形膜（双精度）|
| ktrm6s.f | ? | 6节点三角形膜（单精度）|

**A4. 2D 板弯曲单元（Plate Bending）**
| 文件 | 行数 | 功能 |
|------|:---:|------|
| ktrbsc.f | 657 | 基本弯曲三角形（被 ktrplt/kqdplt 调用）|
| ktrplt.f | 566 | CTRIA 三角形板弯曲 |
| kqdplt.f | 342 | CQUAD 四边形板弯曲 |
| kpltst.f | 48 | 板单元刚度启动 |

**A5. 2D 壳单元（Shell）**
| 文件 | 行数 | 功能 |
|------|:---:|------|
| ktrpld.f | ? | 三角形板壳（双精度）|
| ktrpls.f | ? | 三角形板壳（单精度）|
| ktshld.f | ? | 三角形壳（双精度）|
| ktshls.f | ? | 三角形壳（单精度）|

**A6. 3D 实体单元（Solid）**
| 文件 | 行数 | 功能 |
|------|:---:|------|
| ktetra.f | 384 | CTETRA 四面体实体刚度 |
| ksolid.f | 236 | CHEXA/CWEDGE 实体驱动（分解为四面体）|

**A7. 特殊单元**
| 文件 | 行数 | 功能 |
|------|:---:|------|
| kelas.f | 189 | CELAS 弹簧单元 |
| kconed.f | 904 | 锥壳单元（双精度）|
| kcones.f | ? | 锥壳单元（单精度）|
| kflud2.f | ? | 流体单元（2D）|
| kflud3.f | ? | 流体单元（3D）|
| kflud4.f | ? | 流体单元（4节点）|
| kpanel.f | ? | 面板单元 |
| kslot.f | ? | 槽单元 |
| ktrapr.f | ? | 梯形环单元 |
| ktrirg.f | ? | 三角形环单元 |
| khring.f | ? | 环单元（热）|

**A8. 单元质量矩阵**
| 文件 | 行数 | 功能 |
|------|:---:|------|
| mce1.f - mce1d.f | ~500 | 集中质量计算 |
| mce2.f | ? | 一致质量计算 |
| mflud2-4.f | ? | 流体质量 |
| msolid.f | ? | 实体质量 |
| mring.f | ? | 环单元质量 |
| mtriqd.f | ? | 四边形质量 |
| mtrirg.f | ? | 三角形环质量 |
| mtrplt.f | ? | 板单元质量 |

---

### 类别 B：全局矩阵装配（对应 Zienkiewicz Ch.3）

| 文件 | 行数 | 功能 |
|------|:---:|------|
| sma1.f | ? | 刚度矩阵装配（主驱动）|
| sma1a.f | ? | 刚度矩阵装配（辅助 A）|
| sma1b.f | ? | 刚度矩阵装配（插入子程序）|
| sma2.f | ? | 刚度矩阵装配（第2阶段）|
| sma2a.f | ? | 刚度矩阵装配（辅助）|
| sma3.f | ? | 装配（QUAD4 等参单元）|
| ssg1.f - ssg4.f | ? | 静力载荷生成器 |

---

### 类别 C：方程组求解（对应 Zienkiewicz Ch.3.5-3.6）

**C1. 矩阵分解/求解**
| 文件 | 行数 | 功能 |
|------|:---:|------|
| inverd.f | 124 | 密集矩阵求逆/方程求解（双精度）|
| invers.f | 112 | 密集矩阵求逆/方程求解（单精度）|
| solve.f | 220 | DMAP求解器驱动（调用 sdcomp/fbs）|
| sdcomp.f | 22 | 对称分解入口 |
| decomp.f | 935 | 非对称 LU 分解 |
| fbs.f | 81 | 前后代入（对称）|
| gfbs.f | 357 | 前后代入（一般）|
| invfbs.f | 211 | 逆迭代前后代入 |

**C2. 特征值求解**
| 文件 | 行数 | 功能 |
|------|:---:|------|
| qriter.f | 261 | QR 迭代（三对角矩阵特征值）|
| qriter1.f | 263 | QR 迭代（变体）|
| invpwr.f | 235 | 逆迭代法特征值 |
| invp1.f | 38 | 逆迭代：形成 (K-λM) |
| invp2.f | 81 | 逆迭代：分解 |
| invp3.f | 620 | 逆迭代：迭代主循环 |
| tridi.f | 278 | Givens 三对角化 |
| tridi1.f | 273 | 三对角化（变体）|
| hess1.f | 297 | Hessenberg 特征值变换 |
| reig.f | 353 | 实特征值问题驱动 |
| valvec.f | 136 | 特征值/特征向量驱动 |

**C3. 瞬态分析（时间积分）**
| 文件 | 行数 | 功能 |
|------|:---:|------|
| trd1c.f | 658 | Newmark-β 时间积分主循环 |
| trd1a.f | 149 | 瞬态初始化 |
| trd1d.f | 458 | 瞬态非线性载荷 |
| trd1e.f | 475 | 瞬态输出 |
| timts1.f | 458 | 时间步进（方法1）|
| step.f | 40 | 单步计算 |

---

### 类别 D：后处理/应力恢复（对应 Zienkiewicz Ch.14）

| 文件 | 行数 | 功能 |
|------|:---:|------|
| srod1.f | 118 | 杆单元应力恢复 Phase I |
| srod2.f | 149 | 杆单元应力恢复 Phase II |
| sbar1.f | 549 | 梁单元应力恢复 Phase I |
| sbar2.f | 271 | 梁单元应力恢复 Phase II |
| sdr1.f | 207 | 应力数据恢复驱动 Phase I |
| sdr2.f | 30 | 应力数据恢复驱动 Phase II |
| sdr3.f | 23 | 应力数据恢复驱动 Phase III |

---

### 类别 E：数学工具库（对应 Zienkiewicz Appendix）

| 文件 | 行数 | 功能 |
|------|:---:|------|
| gmmatd.f | ~90 | 通用矩阵乘法（双精度）|
| inverd.f | 124 | 矩阵求逆（双精度）|
| invers.f | 112 | 矩阵求逆（单精度）|
| jacob2.f | 112 | Jacobian 矩阵计算 |
| jacobs.f | 109 | Jacobian（单精度）|
| norm1.f | 29 | 向量归一化 |

---

### 类别 F：坐标变换与几何处理

| 文件 | 行数 | 功能 |
|------|:---:|------|
| transd.f | ? | 坐标系变换（双精度）|
| rotate.f | ? | 旋转矩阵 |
| rotat.f | ? | 旋转 |

---

### 类别 G：I/O 与数据管理（NASTRAN 基础设施）

| 文件 | 行数 | 功能 |
|------|:---:|------|
| read1-7.f | ? | 磁带/文件读取 |
| merge.f | ? | 矩阵合并 |
| partn.f | ? | 矩阵分区 |
| reduce.f | 996 | 矩阵约简 |
| mpyad.f | ? | 矩阵乘加 |

---

## Benchmark 选取与分类映射

基于上述 Zienkiewicz 标准分类，我们的 12 个 benchmark 覆盖情况如下：

| 类别 | Zienkiewicz 章节 | Benchmark | 源文件 |
|------|-----------------|-----------|--------|
| A1. 1D杆 | Ch.5 (Truss) | B1 | krod.f |
| A2. 1D梁 | Ch.5 (Frame) | B2, B9 | kbar.f, kelbow.f |
| A3. 2D膜 | Ch.6 (Plane stress) | B3, B4 | ktrmem.f, kqdmem.f+ktriqd.f |
| A4. 2D板 | Ch.7 (Plate bending) | B5, B6 | ktrplt.f+ktrbsc.f, kqdplt.f+ktrbsc.f |
| A6. 3D实体 | Ch.8 (3D solid) | B7, B8 | ktetra.f, ksolid.f+ktetra.f |
| C1. 方程求解 | Ch.3.5 (Solution) | B10 | inverd.f |
| C2. 特征值 | Ch.12 (Eigenvalue) | B11 | qriter.f |
| D. 后处理 | Ch.14 (Stress recovery) | B12 | srod1.f |

**覆盖说明**：
- 覆盖了有限元**静力学分析**的完整计算链路：单元公式(A) → 求解(C) → 后处理(D)
- 单元类型覆盖了从 1D 到 3D 的全维度梯度
- 装配(B)未作为独立 benchmark：因其深度耦合 NASTRAN 文件系统，且其计算逻辑相对简单（循环调用 SMA1B）
- 瞬态(C3)未纳入：因 trd1c.f 等子程序深度绑定 NASTRAN 外存架构无法独立运行
- 单元类别 A5(壳)/A7(特殊)/A8(质量) 未全部覆盖：留作未来扩展

---

## 统计概要

NASTRAN-95 `mis/` 文件夹共含约 **950 个 .f 文件**，按功能分类约：
- 单元刚度/质量生成 (A)：~150 个文件
- 装配 (B)：~30 个文件
- 求解/特征值/瞬态 (C)：~60 个文件
- 后处理 (D)：~80 个文件
- 数学工具 (E)：~20 个文件
- 坐标变换 (F)：~15 个文件
- I/O 与数据管理 (G)：~200 个文件
- DMAP 执行/控制：~100 个文件
- 绘图/输出格式：~100 个文件
- 其他（空气动力学、热传导、优化等）：~200 个文件

本研究从中选取 12 个代表性子程序作为重构 benchmark，总计 **4,989 行 Fortran 77 源代码**。
