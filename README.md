``` bash
sudo apt install libctl-dev # 必须要
sudo apt install libhdf5-dev # 调试 可选
sudo apt-get -y install libgdsii-dev # 如果需要读取 gds 文件则是必须
```

``` bash
./tests/make_structure  13.5 Air 60 Ru 2.4 $(for ((i=0; i<40; i++)); do echo Si 4 Mo 2.8; done) LTEM 104
```