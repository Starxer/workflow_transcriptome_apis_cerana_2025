#!/bin/bash


while IFS="," read -r old_name new_name; do
    if [ -f "$old_name" ]; then
        mv "$old_name" "$new_name"
        echo "已重命名 $old_name 为 $new_name"
    else
        echo "文件 $old_name 不存在，跳过。"
    fi
done < "$file_list"

