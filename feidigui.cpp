#include <iostream>
#include <vector>
#include <cmath>
#include<stack>
#include <algorithm>

// 节点结构体
struct Node {
    std::vector<double> point;
    Node* left;
    Node* right;

    Node(const std::vector<double>& p) : point(p), left(nullptr), right(nullptr) {}
};

struct Neighbor {
    Node* node;
    double distance;

    bool operator<(const Neighbor& other) const {
        return distance < other.distance;
    }
};
// 计算欧几里得距离
double distance(const std::vector<double>& p1, const std::vector<double>& p2) {
    double dist = 0;
    for (size_t i = 0; i < p1.size(); ++i) {
        dist += std::pow(p1[i] - p2[i], 2);
    }
    return std::sqrt(dist);
}
//构建kd树
Node* buildKdTree(std::vector<std::vector<double>>& points, int depth) {
    if (points.empty()) {
        return nullptr;
    }

    int axis = depth % 3;  // 选择坐标轴

    std::sort(points.begin(), points.end(), [axis](const std::vector<double>& p1, const std::vector<double>& p2) {
        return p1[axis] < p2[axis];
    });

    int medianIndex = points.size() / 2;
    Node* root = new Node(points[medianIndex]);

    std::vector<std::vector<double>> leftPoints(points.begin(), points.begin() + medianIndex);
    std::vector<std::vector<double>> rightPoints(points.begin() + medianIndex + 1, points.end());

    root->left = buildKdTree(leftPoints, depth + 1);
    root->right = buildKdTree(rightPoints, depth + 1);

    return root;
}


// KNN搜索
/*void knnSearch(Node* root, const std::vector<double>& target, int k, std::vector<Node*>& neighbors, int depth) {
    if (root == nullptr) {
        return;
    }

    double dist = distance(root->point, target);

    if (neighbors.size() < k) {
        neighbors.push_back(root);
        std::sort(neighbors.begin(), neighbors.end(), [target](Node* n1, Node* n2) {
            return distance(n1->point, target) < distance(n2->point, target);
        });
    } else if (dist < distance(neighbors.back()->point, target)) {
        neighbors.pop_back();
        neighbors.push_back(root);
        std::sort(neighbors.begin(), neighbors.end(), [target](Node* n1, Node* n2) {
            return distance(n1->point, target) < distance(n2->point, target);
        });
    }

    int axis = depth % 3;
    if (target[axis] < root->point[axis]) {
        knnSearch(root->left, target, k, neighbors, depth + 1);
    } else {
        knnSearch(root->right, target, k, neighbors, depth + 1);
    }

    /*double radius = distance(neighbors.back()->point, target);
    if (std::abs(target[axis] - root->point[axis]) < radius) {
        if (target[axis] < root->point[axis]) {
            knnSearch(root->right, target, k, neighbors, depth + 1);
        } else {
            knnSearch(root->left, target, k, neighbors, depth + 1);
        }
    }
}*/   
/*void knnSearch(Node* root, const std::vector<double>& target, int k, std::vector<Node*>& neighbors, int depth) {
    Node* current = root;
    int currentDepth = depth;

    while (current != nullptr) {
        double dist = distance(current->point, target);

        if (neighbors.size() < k) {
            neighbors.push_back(current);
            std::sort(neighbors.begin(), neighbors.end(), [target](Node* n1, Node* n2) {
                return distance(n1->point, target) < distance(n2->point, target);
            });
        } else if (dist < distance(neighbors.back()->point, target)) {
            neighbors.pop_back();
            neighbors.push_back(current);
            std::sort(neighbors.begin(), neighbors.end(), [target](Node* n1, Node* n2) {
                return distance(n1->point, target) < distance(n2->point, target);
            });
        }

        int axis = currentDepth % 3;
        if (target[axis] < current->point[axis]) {
            current = current->left;
        } else {
            current = current->right;
        }

        currentDepth++;
    }
}*/
void knnSearch(Node* root, const std::vector<double>& target, int k, Neighbor neighbors[], int& neighborCount, int depth) {
    Node* current = root;
    int currentDepth = depth;

    while (current != nullptr) {
        double dist = distance(current->point, target);

        if (neighborCount < k) {
            neighbors[neighborCount++] = {current, dist};
            std::sort(neighbors, neighbors + neighborCount);
        } else if (dist < neighbors[neighborCount - 1].distance) {
            neighbors[neighborCount - 1] = {current, dist};
            std::sort(neighbors, neighbors + neighborCount);
        }

        int axis = currentDepth % 3;
        if (target[axis] < current->point[axis]) {
            current = current->left;
        } else {
            current = current->right;
        }

        currentDepth++;
    }
}



int main() {
    std::vector<std::vector<double>> points = {
        {1, 2, 3},
        {4, 5, 6},
        {5,6,7},
        {7, 8, 9},
        {10, 11, 12},
        {13, 14, 15}
    };

    Node* root = buildKdTree(points, 0);
    Neighbor neighbors[2];
    std::vector<double> target = {5, 6, 7};
    int k = 2;

    int neighborCount = 0;
    knnSearch(root, target, k, neighbors,neighborCount,0);

 std::cout << "最近的 " << k << " 个邻居:" << std::endl;
    for (int i = 0; i < neighborCount; ++i) {
        std::cout << "点: (";
        for (size_t j = 0; j < neighbors[i].node->point.size(); ++j) {
            std::cout << neighbors[i].node->point[j];
            if (j < neighbors[i].node->point.size() - 1) {
                std::cout << ", ";
            }
        }
        std::cout << ") 距离: " << neighbors[i].distance << std::endl;
    }

    return 0;
}
