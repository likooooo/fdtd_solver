``` bash
sudo apt install libctl-dev # 必须要
sudo apt install libhdf5-dev # 调试 可选
sudo apt-get -y install libgdsii-dev # 如果需要读取 gds 文件则是必须
```

``` bash
# lambda = 13.5
# PML thickness = 10
# step = 1
./tests/make_structure  13.5 10 1 Air 60 Ru 2.4 $(for ((i=0; i<40; i++)); do echo Si 4 Mo 2.8; done) LTEM 104
```


``` 问题
0. 网格划分加速方法
1. add_volume_source 填的应该是 shift 相位，而不是 propogate 相位 ？
TODO : 通过改变光源的位置，通过相位进行补偿，得到一个相同的结果
2. 模拟点光源，给一个足够短的持续时间(1 dt) ? 
3. 用 gaussian_src_time 模拟一个点光源，但是 fwidth 定义不对， 应该改用 continuous_src_time 给一个 dt 长度的持续时间
4. structure_chunk::update_condinv 为什么没有找到 conductive?
```

```c++
src_vol_chunkloop : 光源的构建
```


在Meep中，`epsilon` (介电常数 $\epsilon$)、`mu` (磁导率 $\mu$) 和 `conductivity` (电导率 $\sigma$) 都是用来定义**材料电磁特性**的关键参数。它们之间存在密切的关系，共同决定了材料如何与电磁波相互作用，包括光的传播、吸收和反射。

---

### **基本关系**

在时域有限差分 (FDTD) 方法，即 Meep 所基于的方法中，电磁波的传播由麦克斯韦方程组控制。这些方程组在材料中包含以下本构关系：

## 1. 电位移场 (D场) 与电场 (E场)

$$D = \epsilon E$$
$\epsilon$ 是材料的**介电常数**  
它描述了材料在电场作用下储存电能的能力。

## 2. 磁场 (B场) 与磁场强度 (H场)
$$B = \mu H$$
$\mu$ 是材料的**磁导率**  
它描述了材料在磁场作用下储存磁能的能力。


## 3. 传导电流密度 (J场) 与电场 (E场)
$$J = \sigma E$$
$\sigma$ 就是材料的**电导率**  
它描述了材料传导电流的能力，直接导致能量损耗（欧姆损耗）。


#### ** 电导率 ($\sigma$)**

* **直接损耗**：`conductivity` (在 Meep 中通常表示为 `D_conductivity` 用于电场，`B_conductivity` 用于磁场) **直接对应于材料的欧姆损耗**。它在麦克斯韦方程组中表现为电位移场或磁场的虚部。
* **与虚部介电常数/磁导率的关系**：
    * 在频域，电导率 $\sigma_D$ 可以等效为一个虚部的介电常数：
        $\mathrm{Im}(\epsilon) = \frac{\sigma_D}{\omega}$ (其中 $\omega$ 是角频率)
        因此，**`conductivity` 提供了一种直接指定材料损耗的方式**，尤其是在没有复杂频散模型的情况下。
    * 类似的，`B_conductivity` 对应于虚部的磁导率。
* **简化损耗建模**：对于简单的损耗材料，直接设置 `conductivity` 比通过复杂的频散模型来引入 $\mathrm{Im}(\epsilon)$ 更简单直接。尤其是在模拟**金属**时，`conductivity` 是一个非常重要的参数。

---

### **相互作用总结**

* **`epsilon` 和 `mu` 决定了波的传播速度和阻抗**。它们可以具有频率依赖性 (频散) 和复数值 (损耗)。
* **`conductivity` 提供了一种额外的、直接的损耗机制**，它与 `epsilon` 和 `mu` 的虚部作用是等效的，但形式上更直接地与欧姆损耗相关。
* **在 Meep 中，你可以同时设置这些参数**：
    * 你可以为材料指定一个频率无关的 `epsilon` 和 `mu`。
    * 你也可以添加 `susceptibilities` 来定义频率相关的介电响应 (例如 Drude 或 Lorentzian 模型)，这些模型本身就包含了实部和虚部，从而包含损耗。
    * 或者，你可以直接指定 `D_conductivity` 和 `B_conductivity` 来引入欧姆损耗，而无需定义复杂的频散模型。

**通常情况下，对于非磁性介电材料，你主要关注 `epsilon` (及其频散和虚部)。对于金属，`epsilon` 的Drude模型或直接设置 `D_conductivity` 很常见。**

