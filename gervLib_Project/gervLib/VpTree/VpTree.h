#ifndef VPTREE_H
#define VPTREE_H

//#include <queue>
//#include <VpTree/QueueItem.h>
//#include <Dataset.h>
//#include <VpTree/Node/LeafNodeVPTree.h>
//#include <VpTree/Node/DirectorNode.h>
//#include <Pivots.h>

//template <class DType, class DistanceFunction>
//class VpTree{

//    private:
//        Node<DType> *root;
//        Dataset<DType> dataset;
//        Pivot<DType>* algorithm;

//        Node<DType> *buildTree(const bool isCircunscribed,
//                               const double threshold,
//                               const size_t maxElements,
//                               std::vector<BasicArrayObject<DType>> previousPivots,
//                               std::vector<std::vector<double>> distanceVector,
//                               Dataset<DType>* dataset,
//                               DistanceFunction* df);

//        void getBucket(const double threshold,
//                       std::vector<BasicArrayObject<DType>> previousPivots,
//                       std::vector<std::vector<double>> distanceVector,
//                       Dataset<DType> *dataset,
//                       DistanceFunction * df,
//                       Bucket<DType>* bucket);

//        void incrementDistanceVector(const BasicArrayObject<DType> &pivot,
//                                     const BasicArrayObject<DType> &parentPivot,
//                                     std::vector<std::vector<double>> *distanceVector,
//                                     DistanceFunction* df) const;

//        void getNodeHistogram(const Dataset<DType> &dataset,
//                              const BasicArrayObject<DType> &pivot,
//                              DistanceFunction * df,
//                              DirectorNode<DType>* dirNode) const;

//        void balanceTree();

//        LeafNodeVPTree<DType> *generateBalancedLeaf(size_t i,
//                                              size_t j,
//                                              std::vector<std::vector<DirectorNode<DType>*>> heightList);

//        void pushExceededNodes(Node<DType> *node, LeafNodeVPTree<DType> *newLeaf);

//        bool nodeCmpFromPtr(Node<DType> *n1, Node<DType> *n2);
//        auto checkSqPosition(double dist, double mu, double max);
//        void setLeftRightMinMax(const BasicArrayObject<DType> &qElement, DirectorNode<DType> *dirNode, DistanceFunction *df);

//        void kNNIncWalk(const BasicArrayObject<DType> &qElement,
//                      size_t k,
//                      Node<DType> *node,
//                      std::priority_queue<QueueItem<DType>,
//                                          std::vector<QueueItem<DType>>,
//                                          std::greater<QueueItem<DType>>> *queue,
//                      DistanceFunction *df);

//        void kNNWalk(const BasicArrayObject<DType> &qElement,
//                    double &tau,
//                    size_t k,
//                    Node<DType> *node,
//                    std::priority_queue<QueueItem<DType>,
//                                        std::vector<QueueItem<DType>>,
//                                        std::less<QueueItem<DType>>> *queue,
//                    std::vector<std::pair<int, double>> *pivotVec,
//                    DistanceFunction *df);

//    public:
//        VpTree(bool balance,
//               const double paramThreshold,
//               const size_t maxElements,
//               Pivot<DType>* algorithm,
//               Dataset<DType> *dataset,
//               DistanceFunction * df);

//        ~VpTree(){ delete(root); root = nullptr; }

//        void kNN(const BasicArrayObject<DType> &qElement,
//                 size_t k,
//                 Node<DType> *node,
//                 Dataset<DType> *answer,
//                 DistanceFunction * df);

//        void kNNInc(const BasicArrayObject<DType> &qElement,
//                    size_t k,
//                    Node<DType> *node,
//                    Dataset<DType> *answer,
//                    DistanceFunction * df);

//        void rangeQuery(const BasicArrayObject<DType> &qElement,
//                        double radius,
//                        Node<DType> *node,
//                        Dataset<DType> *answer,
//                        DistanceFunction * df);

//        static void splitDataset(const BasicArrayObject<DType> &pivot,
//                          const double mu,
//                          const Dataset<DType> &sourceDataset,
//                          Dataset<DType> *innerDataset,
//                          Dataset<DType> *outerDataset,
//                          DistanceFunction * df);

//        size_t getHeight(Node<DType> *node) const;

//        void getLeafs(Node<DType> *node, std::vector<LeafNodeVPTree<DType>*> *vec) const;

//        void getDirectorNodes(Node<DType> *node, std::vector<DirectorNode<DType> *> *ids);

//        void getDirectorNodesByHeight(Node<DType> *node,
//                                      std::vector<std::vector<DirectorNode<DType>*>> *dirVec,
//                                      size_t height=0) const;

//        Node<DType> *getRoot() const;
//};

////---------------------------------------------------------------------------------------------
////---------------------------------------------------------------------------------------------
////---------------------------------------------------------------------------------------------
////---------------------------------------------------------------------------------------------

//#undef IMPL_TEMPL
//#define IMPL_TEMPL template <class DType, class DistanceFunction>

//IMPL_TEMPL VpTree<DType, DistanceFunction>::VpTree(bool balance,
//                                                   double threshold,
//                                                   size_t maxElements,
//                                                   Pivot<DType>* algorithm,
//                                                   Dataset<DType> *dataset,
//                                                   DistanceFunction *df){

//    this->dataset = *dataset;
//    this->algorithm = algorithm;

//    std::vector<BasicArrayObject<DType>> previousPivots;
//    std::vector<std::vector<double>> distanceVector;

//    root = buildTree(false, threshold, maxElements, previousPivots, distanceVector, dataset, df/*, progressBar*/);

//    if (balance) { balanceTree(); }
//}

//IMPL_TEMPL Node<DType> *VpTree<DType, DistanceFunction>::buildTree(const bool isCircunscribed,
//                                                                   const double threshold,
//                                                                   const size_t maxElements,
//                                                                   std::vector<BasicArrayObject<DType>> previousPivots,
//                                                                   std::vector<std::vector<double>> distanceVector,
//                                                                   Dataset<DType> *dataset,
//                                                                   DistanceFunction *df){

////    const DType pivot = PivotSelection::getPivot(*dataset, algorithm, df);
//    this->algorithm->generatePivots(dataset, df, 1);
//    const BasicArrayObject<DType> pivot = *this->algorithm->getPivot(0);

//    const double mu = dataset->medianDistance(pivot, df);

//    if (mu > threshold && dataset->getSize() > maxElements){

//        DirectorNode<DType> *dirNode = new DirectorNode<DType>;

//        dirNode->setPivot(pivot);
//        dirNode->setMu(mu);
//        dirNode->setCoverage(Dataset<DType>::getMax(*dataset, pivot, df));

//        incrementDistanceVector(pivot, previousPivots.back(), &distanceVector, df);
//        previousPivots.push_back(pivot);

//        Dataset<DType> *innerDataset = new Dataset<DType>;
//        Dataset<DType> *outerDataset = new Dataset<DType>;
//        splitDataset(pivot, mu, *dataset, innerDataset, outerDataset, df);

//        delete(dataset);

//        if (innerDataset->getSize() > 0){
//            dirNode->setLeftNode(buildTree(true, threshold, maxElements, previousPivots, distanceVector,
//                                           innerDataset, df/*, progressBar*/));
//        } else {
//            delete(innerDataset);
//            dirNode->setLeftNode(nullptr);
//        }

//        if (outerDataset->getSize() > 0){
//            dirNode->setRightNode(buildTree(false, threshold, maxElements, previousPivots, distanceVector,
//                                            outerDataset, df/*, progressBar*/));
//        } else {
//            delete(outerDataset);
//            dirNode->setRightNode(nullptr);
//        }

//        return dirNode;

//    } else {

//        LeafNodeVPTree<BasicArrayObject<DType>> *leaf = new LeafNodeVPTree<BasicArrayObject<DType>>;

//        leaf->setCircunscribed(isCircunscribed);
//        leaf->setPivot(pivot);
//        leaf->setMu(mu);

//        leaf->setCoverage(Dataset<DType>::getMax(*dataset, pivot, df));

//        incrementDistanceVector(pivot, previousPivots.back(), &distanceVector, df);
//        leaf->setDistanceVector(distanceVector);

//        previousPivots.push_back(pivot);
//        leaf->setPreviousPivots(previousPivots);

//        for (size_t i = 0; i < dataset->getSize(); i++){

//            leaf->push_back(Pair<BasicArrayObject<DType>>(dataset->operator [](i), distanceVector.back()));

//        }

//        delete(dataset);

//        return leaf;
//    }
//}

//IMPL_TEMPL void VpTree<DType, DistanceFunction>::getBucket(const double threshold,
//                                                           std::vector<BasicArrayObject<DType>> previousPivots,
//                                                           std::vector<std::vector<double>> distanceVector,
//                                                           Dataset<DType> *dataset,
//                                                           DistanceFunction *df,
//                                                           Bucket<DType> *bucket){

//    srand(time(nullptr));
//    const BasicArrayObject<DType> pivot = dataset->operator [](rand()% dataset->getCardinality());

//    const double mu = dataset->medianDistance(pivot, df);

//    if ((mu > threshold) && (dataset->getCardinality() > 1)){

//        // Split the dataset in innerDataset (< mu) and outerDataset (>= mu).
//        Dataset<DType> *innerDataset = new Dataset<DType>;
//        Dataset<DType> *outerDataset = new Dataset<DType>;
//        splitDataset(pivot, mu, *dataset, innerDataset, outerDataset, df);

//        // Deletes the dataset to save memory.
//        delete(dataset);

//        // Increments the DistanceVector with a vector of distances for this node.
//        incrementDistanceVector(pivot, previousPivots.back(), &distanceVector, df);

//        bucket->push_back(Pair<DType>(pivot, distanceVector.back()));

//        previousPivots.push_back(pivot);

//        // Recursively calls itself with the inner part of the dataset.
//        if (innerDataset->getSize() > 0){
//            getBucket(threshold, previousPivots, distanceVector, innerDataset, df, bucket/*, progressBar*/);
//        } else {
//            delete(innerDataset);
//        }

//        // Recursively calls itself with the outer part of the dataset.
//        if (outerDataset->getSize() > 0){
//            getBucket(threshold, previousPivots, distanceVector, outerDataset, df, bucket/*, progressBar*/);
//        } else {
//            delete(outerDataset);
//        }

//    } else {
//        delete(dataset);

//        incrementDistanceVector(pivot, previousPivots.back(), &distanceVector, df);

//        bucket->push_back(Pair<DType>(pivot, distanceVector.back()));
//    }

//}

//IMPL_TEMPL void VpTree<DType, DistanceFunction>::splitDataset(const BasicArrayObject<DType> &pivot,
//                                                              const double mu,
//                                                              const Dataset<DType> &sourceDataset,
//                                                              Dataset<DType> *innerDataset,
//                                                              Dataset<DType> *outerDataset,
//                                                              DistanceFunction *df){

//    const size_t cardinality = sourceDataset.getCardinality();

////    const uint32_t pivotOID = pivot.getOID();

//    for (size_t i = 0; i < cardinality; i++){

//        const BasicArrayObject<DType> &fv = sourceDataset.getFeatureVector(i);

////        if (pivotOID != fv.getOID()){
//            if (df->getDistance(pivot, fv) < mu)
//                innerDataset->push_back(fv);
//            else
//                outerDataset->push_back(fv);
////        }

//    }

//    innerDataset->setCardinality(innerDataset->getSize());
//    outerDataset->setCardinality(outerDataset->getSize());

//    innerDataset->setDimensionality(sourceDataset.getDimensionality());
//    outerDataset->setDimensionality(sourceDataset.getDimensionality());
//}

//IMPL_TEMPL void VpTree<DType, DistanceFunction>::incrementDistanceVector(const BasicArrayObject<DType> &pivot, const BasicArrayObject<DType> &parentPivot,
//                                                                         std::vector<std::vector<double>> *distanceVector,
//                                                                         DistanceFunction *df) const{

//    std::vector<double> distances;

//    if (distanceVector->begin() != distanceVector->end()){

//        distances = distanceVector->back();

//        const double distanceToParent = df->getDistance(pivot, parentPivot);

//        // Increments every element of distances by distanceToParent.
//        std::for_each(distances.begin(), distances.end(), [distanceToParent](double &dist){
//            dist += distanceToParent;
//        });

//        distances.push_back(distanceToParent);
//    }
//    distanceVector->push_back(distances);

//}

//IMPL_TEMPL size_t VpTree<DType, DistanceFunction>::getHeight(Node<DType> *node) const{

//    if (node == nullptr) return -1;

//    DirectorNode<DType> *dirNode = dynamic_cast<DirectorNode<DType>*>(node);

//    if (dirNode != nullptr)
//        return 1 + std::max(getHeight(dirNode->getLeftNode()), getHeight(dirNode->getRightNode()));
//    else
//        return 1;
//}

//IMPL_TEMPL void VpTree<DType, DistanceFunction>::getLeafs(Node<DType> *node, std::vector<LeafNodeVPTree<DType>*>*leafVec) const{

//    if (node != nullptr){

//        const DirectorNode<DType> *dirNode = dynamic_cast<DirectorNode<DType>*>(node);

//        if (dirNode != nullptr){
//            getLeafs(dirNode->getLeftNode(), leafVec);
//            getLeafs(dirNode->getRightNode(), leafVec);
//        } else {
//            leafVec->push_back(dynamic_cast<LeafNodeVPTree<DType>*>(node));
//        }

//    }

//}

//IMPL_TEMPL void VpTree<DType, DistanceFunction>::getDirectorNodes(Node<DType> *node, std::vector<DirectorNode<DType>*> *ids){

//    if (node != nullptr){

//        DirectorNode<DType> *dirNode = dynamic_cast<DirectorNode<DType>*>(node);

//        if (dirNode != nullptr){

//            ids->push_back(dirNode);

//            getDirectorNodes(dirNode->getLeftNode(), ids);
//            getDirectorNodes(dirNode->getRightNode(), ids);
//        }

//    }

//}

//IMPL_TEMPL void VpTree<DType, DistanceFunction>::getDirectorNodesByHeight(Node<DType> *node,
//        std::vector<std::vector<DirectorNode<DType>*>> *dirVec,
//        size_t height) const{

//    if (node != nullptr){

//        DirectorNode<DType> *dirNode = dynamic_cast<DirectorNode<DType>*>(node);

//        if (dirNode != nullptr){

//            dirVec->at(height).push_back(dirNode);

//            getDirectorNodesByHeight(dirNode->getLeftNode(), dirVec, height+1);
//            getDirectorNodesByHeight(dirNode->getRightNode(), dirVec, height+1);
//        }

//    }
//}

//IMPL_TEMPL LeafNodeVPTree<DType> *VpTree<DType, DistanceFunction>::generateBalancedLeaf(size_t i, size_t j, std::vector<std::vector<DirectorNode<DType>*>> heightList){

//    LeafNodeVPTree<DType> *newLeaf = new LeafNodeVPTree<DType>;

//    newLeaf->setPivot(heightList[i+1][j]->getPivot());
//    newLeaf->setMu(heightList[i+1][j]->getMu());
//    newLeaf->setCircunscribed(true);

//    pushExceededNodes(heightList[i+1][j]->getLeftNode(), newLeaf);
//    pushExceededNodes(heightList[i+1][j]->getRightNode(), newLeaf);

//    return newLeaf;
//}

//IMPL_TEMPL void VpTree<DType, DistanceFunction>::balanceTree(){

//    int height = getHeight(root);
//    std::vector<std::vector<DirectorNode<DType>*>> heightList(height);
//    getDirectorNodesByHeight(root, &heightList);

//    for (size_t i = 0; i < heightList.size()-1; i++){

//        // If next height is less than double the size of actual height, the next height is unbalanced.
//        if (heightList[i].size()*2 > heightList[i+1].size()){

//            // For every node in the unbalanced height.
//            for (size_t j = 0; j < heightList[i+1].size(); j++){

//                // Search for nodes pointing to the unbalanced node and points it to the new node.
//                for (size_t k = 0; k < heightList[i].size(); k++){

//                    if (heightList[i][k]->getLeftNode() == heightList[i+1][j])
//                        heightList[i][k]->setLeftNode(generateBalancedLeaf(i, j, heightList));

//                    if (heightList[i][k]->getRightNode() == heightList[i+1][j])
//                        heightList[i][k]->setRightNode(generateBalancedLeaf(i, j, heightList));
//                }
//            }

//            break;
//        }
//    }
//}

//IMPL_TEMPL void VpTree<DType, DistanceFunction>::pushExceededNodes(Node<DType> *node,
//                                                                   LeafNodeVPTree<DType> *newLeaf){

//    if (node == nullptr) return;

//    newLeaf->push_exceeded_node(node);

//    DirectorNode<DType> *dirNode = dynamic_cast<DirectorNode<DType>*>(node);

//    if (dirNode != nullptr){

//        pushExceededNodes(dirNode->getLeftNode(), newLeaf);
//        pushExceededNodes(dirNode->getRightNode(), newLeaf);

//    } else {

//        LeafNodeVPTree<DType> *leafNode = static_cast<LeafNodeVPTree<DType>*>(node);

//        for (size_t i = 0; i < leafNode->getBucket().numberOfElements(); i++){
//            newLeaf->push_back(leafNode->getPair(i));
//        }

//    }
//}

///**
// *@brief A getter to the root.
// *@returns a pointer to the root of the tree
// */
//IMPL_TEMPL Node<DType> *VpTree<DType, DistanceFunction>::getRoot() const{ return root; }

//IMPL_TEMPL void VpTree<DType, DistanceFunction>::rangeQuery(const BasicArrayObject<DType> &qElement,
//                                                            double radius,
//                                                            Node<DType> *node,
//                                                            Dataset<DType> *answer,
//                                                            DistanceFunction *df){


//     if (node == nullptr) return;
//     node->wasVisited = true;

//     if (df->getDistance(qElement, node->getPivot()) <= radius)
//         answer->push_back(node->getPivot());

//     LeafNodeVPTree<DType> *leaf = dynamic_cast<LeafNodeVPTree<DType>*>(node);

//     if (leaf != nullptr) {

//         for (size_t i = 0; i < leaf->numberOfElements(); i++){

//             if (df->getDistance(qElement, leaf->getPair(i).first()) <= radius)
//                 answer->push_back(leaf->getPair(i).first());
//         }

//     } else {

//        double dist = df->getDistance(node->getPivot(), qElement);

//        DirectorNode<DType> *dirNode = dynamic_cast<DirectorNode<DType>*>(node);

//        if (dist < node->getMu()) {

//            if ((dist - radius) <= node->getMu())
//                rangeQuery(qElement, radius, dirNode->getLeftNode(), answer, df);

//            if ((dist + radius) >= node->getMu())
//                rangeQuery(qElement, radius, dirNode->getRightNode(), answer, df);

//        } else {

//            if ((dist + radius) >= node->getMu())
//                rangeQuery(qElement, radius, dirNode->getRightNode(), answer, df);

//            if ((dist - radius) <= node->getMu())
//                rangeQuery(qElement, radius, dirNode->getLeftNode(), answer, df);
//        }

//     }

//}

//enum PositionRelativeToPartition{INSIDE, RING, OUTSIDE};

//IMPL_TEMPL auto VpTree<DType, DistanceFunction>::checkSqPosition(double dist, double mu, double max) {
//    if (dist < mu)
//        return PositionRelativeToPartition::INSIDE;
//    else if (dist <= max)
//        return PositionRelativeToPartition::RING;
//    else
//        return PositionRelativeToPartition::OUTSIDE;
//}

//IMPL_TEMPL void VpTree<DType, DistanceFunction>::setLeftRightMinMax(const BasicArrayObject<DType> &qElement,
//                                                                    DirectorNode<DType> *dirNode,
//                                                                    DistanceFunction *df){


//    double dist =  df->getDistance(qElement, dirNode->getPivot());

//    PositionRelativeToPartition sqPosition = checkSqPosition(dist, dirNode->getMu(), dirNode->getCoverage());

//    // sq dentro bola
//    if (sqPosition == PositionRelativeToPartition::INSIDE) {

//        // Min do Left
//        dirNode->getLeftNode()->setDMin(0.0);

//        // Max do Left
//        if (dirNode->getDMax() > dirNode->getMu() + dist)
//            dirNode->getLeftNode()->setDMax(dirNode->getMu() + dist);
//        else
//            dirNode->getLeftNode()->setDMax(dirNode->getDMax());

//        // Min do Right
//        dirNode->getRightNode()->setDMin(dirNode->getMu() - dist);

//        // Max do Right
//        if (dirNode->getDMax() > dirNode->getCoverage() + dist)
//            dirNode->getRightNode()->setDMax(dirNode->getCoverage() + dist);
//        else
//            dirNode->getRightNode()->setDMax(dirNode->getDMax());
//    }
//    // sq no anel
//    else if (sqPosition == PositionRelativeToPartition::RING) {

//        // Min do Left
//        dirNode->getLeftNode()->setDMin(dist - dirNode->getMu());

//        // Max do Left
//        if (dirNode->getDMax() > dist + dirNode->getMu())
//            dirNode->getLeftNode()->setDMax(dist + dirNode->getMu());
//        else
//            dirNode->getLeftNode()->setDMax(dirNode->getDMax());

//        // Min do Right
//        dirNode->getRightNode()->setDMin(0.0);

//        // Max do Right
//        if (dirNode->getDMax() > dist + dirNode->getCoverage())
//            dirNode->getRightNode()->setDMax(dist + dirNode->getCoverage());
//        else
//            dirNode->getRightNode()->setDMax(dirNode->getDMax());
//    }
//    // sq fora da partição
//    else {

//        // Min do Left
//        dirNode->getLeftNode()->setDMin(dist - dirNode->getMu());

//        // Max do Left
//        if (dirNode->getDMax() > dist + dirNode->getMu())
//            dirNode->getLeftNode()->setDMax(dist + dirNode->getMu());
//        else
//            dirNode->getLeftNode()->setDMax(dirNode->getDMax());

//        // Min do Right
//        dirNode->getRightNode()->setDMin(dist - dirNode->getCoverage());

//        // Max do Right
//        if (dirNode->getDMax() > dist + dirNode->getCoverage())
//            dirNode->getRightNode()->setDMax(dist + dirNode->getCoverage());
//        else
//            dirNode->getRightNode()->setDMax(dirNode->getDMax());
//    }
//}

//template <class DType>
//class NodeCmpFromPtr{

//public:
//    bool operator() (Node<DType> *node1, Node<DType> *node2){
//        return node1->operator>(*node2);
//    }
//};



//IMPL_TEMPL void VpTree<DType, DistanceFunction>::kNNIncWalk(const BasicArrayObject<DType> &qElement,
//                                                            size_t k,
//                                                            Node<DType> *tree,
//                                                            std::priority_queue<QueueItem<DType>,
//                                                                                std::vector<QueueItem<DType>>,
//                                                                                std::greater<QueueItem<DType>>> *resultQueue,
//                                                            DistanceFunction *df){

//    std::priority_queue<Node<DType> *,
//                        std::vector<Node<DType> *>,
//                        NodeCmpFromPtr<DType>> nodeQueue;

//    std::priority_queue<QueueItem<DType>,
//                        std::vector<QueueItem<DType>>,
//                        std::greater<QueueItem<DType>> > elementQueue;

//    tree->setDMin(0.0);
//    tree->setDMax(std::numeric_limits<double>::infinity());

//    nodeQueue.push(tree);

//    // Não completou resultado && há nós ou elementos a serem obtidos.
//    while (resultQueue->size() < k && !(nodeQueue.empty() && elementQueue.empty())){

//        // Se entrar aqui, é garantido que o nodeQueue tem nós e que o elementQueue está vazio.
//        if (elementQueue.empty()) {

//            // Popa a partição
//            Node<DType> *node = nodeQueue.top();
//            nodeQueue.pop();

//            // Se o cast falhar, nó não é do tipo Leaf e retorna nullptr.
//            LeafNodeVPTree<DType> *leaf = dynamic_cast<LeafNodeVPTree<DType> *>(node);
//            if (leaf != nullptr) {

//                // Insere os elementos do folha no elementQueue
//                for (int i = 0; i < leaf->numberOfElements(); i++){
//                    BasicArrayObject<DType> leafElement = leaf->getPair(i).first();
//                    double dist = df->getDistance(leafElement, qElement);
//                    QueueItem<DType> qi = QueueItem<DType>(dist, leafElement);
//                    elementQueue.push(qi);
//                }

//                // Seta variável para verificarmos quais partições foram visitadas.
//                node->wasVisited = true;

//            } else {

//                DirectorNode<DType> *dirNode = dynamic_cast<DirectorNode<DType> *>(node);

//                // Seta min/max dos filhos
//                setLeftRightMinMax(qElement, dirNode, df);

//                // Insere os filhos no nodeQueue
//                nodeQueue.push(dirNode->getLeftNode());
//                nodeQueue.push(dirNode->getRightNode());

//            }

//            // Se entrar aqui, é garantido que o elementQueue tem elementos e que o nodeQueue tem nós,
//            // assim como que o dMin do topo do nodeQueue é menor do que a distância no topo do elementQueue.
//        } else if (!nodeQueue.empty() &&
//                   nodeQueue.top()->getDMin() < elementQueue.top().dist) {

//            // Popa a partição
//            Node<DType> *node = nodeQueue.top();
//            nodeQueue.pop();

//            // Checa se é folha
//            // Se o cast falhar, nó não é do tipo Leaf e retorna nullptr.
//            LeafNodeVPTree<DType> *leaf = dynamic_cast<LeafNodeVPTree<DType>*>(node);
//            if (leaf != nullptr) {

//                // Insere os elementos do folha no elementQueue
//                for (size_t i = 0; i < leaf->numberOfElements(); i++){
//                    BasicArrayObject<DType> leafElement = leaf->getPair(i).first();
//                    double dist = df->getDistance(leafElement, qElement);
//                    QueueItem<DType> qi = QueueItem<DType>(dist, leafElement);
//                    elementQueue.push(qi);
//                }

//                // Seta variável para verificarmos quais partições foram visitadas.
//                node->wasVisited = true;

//            } else {

//                DirectorNode<DType> *dirNode = dynamic_cast<DirectorNode<DType>*>(node);

//                // Seta min/max dos filhos
//                setLeftRightMinMax(qElement, dirNode, df);

//                // Insere os filhos no nodeQueue
//                nodeQueue.push(dirNode->getLeftNode());
//                nodeQueue.push(dirNode->getRightNode());

//            }

//            // Se entrar aqui, é garantido ou que __nodeQueue está vazio, mas o elementQueue não__,
//            // ou que __a distância do elemento no topo do elementQueue é menor do que o dMin do topo do nodeQueue__
//        } else {
//            resultQueue->push(elementQueue.top());
//            elementQueue.pop();
//        }
//    }
//}

//IMPL_TEMPL void VpTree<DType, DistanceFunction>::kNNInc(const BasicArrayObject<DType> &qElement,
//                                                     size_t k,
//                                                     Node<DType> *node,
//                                                     Dataset<DType> *answer,
//                                                     DistanceFunction *df){

//    std::priority_queue<QueueItem<DType>,
//                        std::vector<QueueItem<DType>>,
//                        std::greater<QueueItem<DType>>> queue;

//    kNNIncWalk(qElement, k, node, &queue, df);

////    std::cout << "(kNNInc) resultQueue.size(): " << queue.size() << std::endl;
//    answer->reserve(k);
//    for (size_t i = 0; i < k; i++){
//        answer->push_back(queue.top().featureVector);
////        std::cout << "(kNNInc) Distancia do " << i << " vizinho mais proximo:" << queue.top().dist << std::endl;
//        queue.pop();
//    }
//}


//IMPL_TEMPL void VpTree<DType, DistanceFunction>::kNNWalk(const BasicArrayObject<DType> &qElement,
//                                                         double &tau,
//                                                         size_t k,
//                                                         Node<DType> *node,
//                                                         std::priority_queue<QueueItem<DType>,
//                                                                             std::vector<QueueItem<DType>>,
//                                                                             std::less<QueueItem<DType>>> *queue,
//                                                         std::vector<std::pair<int, double>> *pivotVec,
//                                                         DistanceFunction *df){

//    node->wasVisited = true;

//    LeafNodeVPTree<DType> *leaf = dynamic_cast<LeafNodeVPTree<DType>*>(node);

//    if (leaf != nullptr) {

//        for (size_t i = 0; i < leaf->numberOfElements(); i++){

//            BasicArrayObject<DType> leafElement = leaf->getPair(i).first();
//            double dist = 0.0;

//            bool found = 0;
//            for (size_t j = 0; j < pivotVec->size(); j++){
//                if (pivotVec->operator[](j).first == leafElement.getOID()){
//                    dist = pivotVec->operator[](j).second;
//                    found = 1;
//                    break;
//                }
//            }

//            if (!found)
//                dist = df->getDistance(qElement, leafElement);

//            queue->push(QueueItem<DType>(dist, leafElement));

//            if (queue->size() > k)
//                queue->pop();
//        }

//    } else {

//        DirectorNode<DType> *dirNode = dynamic_cast<DirectorNode<DType>*>(node);

//        double dist = df->getDistance(dirNode->getPivot(), qElement);

//        pivotVec->push_back(std::pair<int, double>(dirNode->getPivot().getOID(), dist));

////        if (dist <= tau){

////            queue->push(QueueItem<DType>(dist, dirNode->getPivot()));

////            while (queue->size() > k)
////                queue->pop();

////            if (queue->size() == k)
////                tau = queue->top().dist;
////        }

//        if (dist < node->getMu()){

//            if ((dist - tau <= dirNode->getMu()) && dirNode->getLeftNode() != nullptr)
//                kNNWalk(qElement, tau, k, dirNode->getLeftNode(), queue, pivotVec, df);

//            if ((dist + tau >= dirNode->getMu()) && dirNode->getRightNode() != nullptr)
//                kNNWalk(qElement, tau, k, dirNode->getRightNode(), queue, pivotVec, df);

//        } else {

//            if ((dist + tau >= dirNode->getMu()) && dirNode->getRightNode() != nullptr)
//                kNNWalk(qElement, tau, k, dirNode->getRightNode(), queue, pivotVec, df);

//            if ((dist - tau <= node->getMu()) && dirNode->getLeftNode() != nullptr)
//                kNNWalk(qElement, tau, k, dirNode->getLeftNode(), queue, pivotVec, df);

//        }

//    }
//}

//IMPL_TEMPL void VpTree<DType, DistanceFunction>::kNN(const BasicArrayObject<DType> &qElement,
//                                                     size_t k,
//                                                     Node<DType> *node,
//                                                     Dataset<DType> *answer,
//                                                     DistanceFunction *df){

//    std::priority_queue<QueueItem<DType>,
//                        std::vector<QueueItem<DType>>,
//                        std::less<QueueItem<DType>>> queue;

//    double tau = std::numeric_limits<double>::infinity();
//    kNNWalk(qElement, tau, k, node, &queue, df);


//    std::priority_queue<QueueItem<DType>,
//                        std::vector<QueueItem<DType>>,
//                        std::greater<QueueItem<DType>>> reverse;

//    for (size_t i = 0; i < k; i++){
//        reverse.push(queue.top());
//        queue.pop();
//    }

//    //std::cout << "(kNN Padrão) resultQueue.size(): " << reverse.size() << std::endl;

//    answer->reserve(k);
//    for (size_t i = 0; i < k; i++){
//        answer->push_back(reverse.top().featureVector);
//        //std::cout << "(kNN Padrao) Distancia do " << i << " vizinho mais proximo:" << reverse.top().dist << std::endl;
//        reverse.pop();
//    }


//}

////IMPL_TEMPL void VpTree<DType, DistanceFunction>::updateProgress(QProgressBar *progressBar){
//////    ++progress;
//////    if (progressBar->maximum() < SHRT_MAX)
//////        emit progressBarValueChanged(progress);
//////    else if (progress % SHRT_MAX == 0)
//////        emit progressBarValueChanged(progress/SHRT_MAX);
////}

#include <queue>
#include <VpTree/QueueItem.h>
#include <Dataset.h>
#include <VpTree/Node/LeafNodeVPTree.h>
#include <VpTree/Node/DirectorNode.h>
#include <Pivots.h>
#include <config_spb.h>
#include <MemoryManagerUtils.h>

template <class type, class DistanceFunction>
class VpTree{

    private:
        Node<BasicArrayObject<type>> *root;
        Dataset<type> dataset;
        Pivot<type>* algorithm;
        size_t progress = 0;
        size_t leafNodeAccess = 0;

        Node<BasicArrayObject<type>> *buildTree(const bool isCircunscribed,
                               const double threshold,
                               const size_t maxElements,
                               std::vector<BasicArrayObject<type>> previousPivots,
                               std::vector<std::vector<double>> distanceVector,
                               Dataset<type> *dataset,
                               DistanceFunction * df/*,
                               QProgressBar *progressBar*/);

        void getBucket(const double threshold,
                       std::vector<BasicArrayObject<type>> previousPivots,
                       std::vector<std::vector<double>> distanceVector,
                       Dataset<type> *dataset,
                       DistanceFunction * df,
                       Bucket<BasicArrayObject<type>>* bucket/*,
                       QProgressBar *progressBar*/);

        void incrementDistanceVector(const BasicArrayObject<type> &pivot,
                                     const BasicArrayObject<type> &parentPivot,
                                     std::vector<std::vector<double>> *distanceVector,
                                     DistanceFunction * df) const;

//        void getNodeHistogram(const Dataset<type> &dataset,
//                              const BasicArrayObject<type> &pivot,
//                              DistanceFunction * df,
//                              DirectorNode<BasicArrayObject<type>>* dirNode) const;

        void balanceTree();

        LeafNodeVPTree<BasicArrayObject<type>> *generateBalancedLeaf(size_t i,
                                              size_t j,
                                              std::vector<std::vector<DirectorNode<BasicArrayObject<type>>*>> heightList);

        void pushExceededNodes(Node<BasicArrayObject<type>> *node, LeafNodeVPTree<BasicArrayObject<type>> *newLeaf);

        bool nodeCmpFromPtr(Node<BasicArrayObject<type>> *n1, Node<BasicArrayObject<type>> *n2);
        auto checkSqPosition(double dist, double mu, double max);
        void setLeftRightMinMax(const BasicArrayObject<type> &qElement, DirectorNode<BasicArrayObject<type>> *dirNode, DistanceFunction *df);

        void kNNIncWalk(const BasicArrayObject<type> &qElement,
                      size_t k,
                      Node<BasicArrayObject<type>> *node,
                      std::priority_queue<QueueItem<BasicArrayObject<type>>,
                                          std::vector<QueueItem<BasicArrayObject<type>>>,
                                          std::greater<QueueItem<BasicArrayObject<type>>>> *queue,
                      DistanceFunction *df);

        void kNNWalk(const BasicArrayObject<type> &qElement,
                    double &tau,
                    size_t k,
                    Node<BasicArrayObject<type>> *node,
                    std::priority_queue<QueueItem<BasicArrayObject<type>>,
                                        std::vector<QueueItem<BasicArrayObject<type>>>,
                                        std::less<QueueItem<BasicArrayObject<type>>>> *queue,
                    std::vector<std::pair<int, double>> *pivotVec,
                    DistanceFunction *df);

        void initDisk()
        {

            baseFilePath = "../VpTree/vp_files";
            setBaseFilePath("VPfiles");

            std::queue<Node<BasicArrayObject<type>>*> queue;
            queue.push(getRoot());
            size_t id = 0;

            while(!queue.empty())
            {

                Node<BasicArrayObject<type>>* curr = queue.front();
                queue.pop();

                LeafNodeVPTree<BasicArrayObject<type>> *leaf = dynamic_cast<LeafNodeVPTree<BasicArrayObject<type>> *>(curr);

                if(leaf != nullptr)
                {

                    leaf->setPageID(id++);
                    saveLeafNode(curr);

                }
                else
                {

                    DirectorNode<BasicArrayObject<type>> *dirNode = dynamic_cast<DirectorNode<BasicArrayObject<type>> *>(curr);
                    if(dirNode->getLeftNode() != nullptr) queue.push(dirNode->getLeftNode());
                    if(dirNode->getRightNode() != nullptr) queue.push(dirNode->getRightNode());

                }

            }

        }

        void saveLeafNode(Node<BasicArrayObject<type>>* curr)
        {

            std::vector<BasicArrayObject<type>> data;
            LeafNodeVPTree<BasicArrayObject<type>> *leaf = dynamic_cast<LeafNodeVPTree<BasicArrayObject<type>> *>(curr);

            for (int i = 0; i < leaf->numberOfElements(); i++)
            {

                data.push_back(leaf->getPair(i).first());

            }

            leaf->resetBucket();

            Dataset<type>* datasetLeaf = new Dataset<type>(data, data.size(), data[0].size());
            write_dataset_to_disk(datasetLeaf, curr->getPageID());
            delete datasetLeaf;

        }

        Dataset<type>* readLeafNode(Node<BasicArrayObject<type>>* curr)
        {

            Dataset<type>* datasetLeaf = new Dataset<type>();
            read_dataset_from_disk(datasetLeaf, curr->getPageID());
            return datasetLeaf;

        }

//        void updateProgress(QProgressBar *progressBar);

    public:
        VpTree(bool balance,
               const double paramThreshold,
               const size_t maxElements,
               Pivot<type>* algorithm,
               Dataset<type> *dataset,
               DistanceFunction * df/*,
               QProgressBar *progressBar*/);

        ~VpTree(){ delete(root); root = nullptr; }

        size_t getLeafNodeAccess(){ return leafNodeAccess; }

        void kNN(const BasicArrayObject<type> &qElement,
                 size_t k,
                 Node<BasicArrayObject<type>> *node,
                 Dataset<type> *answer,
                 DistanceFunction * df);

        void kNNInc(const BasicArrayObject<type> &qElement,
                    size_t k,
                    Node<BasicArrayObject<type>> *node,
                    Dataset<type> *answer,
                    DistanceFunction * df);

        void rangeQuery(const BasicArrayObject<type> &qElement,
                        double radius,
                        Node<BasicArrayObject<type>> *node,
                        Dataset<type> *answer,
                        DistanceFunction * df);

        static void splitDataset(const BasicArrayObject<type> &pivot,
                          const double mu,
                          const Dataset<type> &sourceDataset,
                          Dataset<type> *innerDataset,
                          Dataset<type> *outerDataset,
                          DistanceFunction * df);

        size_t getHeight(Node<BasicArrayObject<type>> *node) const;

        void getLeafs(Node<BasicArrayObject<type>> *node, std::vector<LeafNodeVPTree<BasicArrayObject<type>>*> *vec) const;

        void getDirectorNodes(Node<BasicArrayObject<type>> *node, std::vector<DirectorNode<BasicArrayObject<type>> *> *ids);

        void getDirectorNodesByHeight(Node<BasicArrayObject<type>> *node,
                                      std::vector<std::vector<DirectorNode<BasicArrayObject<type>>*>> *dirVec,
                                      size_t height=0) const;

        Node<BasicArrayObject<type>> *getRoot() const;
};

#undef IMPL_TEMPL
#define IMPL_TEMPL template <class type, class DistanceFunction>

IMPL_TEMPL VpTree<type, DistanceFunction>::VpTree(bool balance,
                                                   double threshold,
                                                   size_t maxElements,
                                                   Pivot<type>* algorithm,
                                                   Dataset<type> *dataset,
                                                   DistanceFunction *df/*,
                                                   QProgressBar *progressBar*/){

    this->dataset = *dataset;
    this->algorithm = algorithm;

    //connect(this, SIGNAL(progressBarValueChanged(int)), progressBar, SLOT(setValue(int)));

    std::vector<BasicArrayObject<type>> previousPivots;
    std::vector<std::vector<double>> distanceVector;

    root = buildTree(false, threshold, maxElements, previousPivots, distanceVector, dataset, df/*, progressBar*/);

    if (balance) { balanceTree(); }

    initDisk();

}

IMPL_TEMPL Node<BasicArrayObject<type>> *VpTree<type, DistanceFunction>::buildTree(const bool isCircunscribed,
                                                                   const double threshold,
                                                                   const size_t maxElements,
                                                                   std::vector<BasicArrayObject<type>> previousPivots,
                                                                   std::vector<std::vector<double>> distanceVector,
                                                                   Dataset<type> *dataset,
                                                                   DistanceFunction *df/*,
                                                                   QProgressBar *progressBar*/){

    //updateProgress(progressBar);

    this->algorithm->generatePivots(dataset, df, 1);
    const BasicArrayObject<type> pivot = *algorithm->getPivot(0);

    const double mu = dataset->medianDistance(pivot, df);

    if (mu > threshold && dataset->getSize() > maxElements){

        DirectorNode<BasicArrayObject<type>> *dirNode = new DirectorNode<BasicArrayObject<type>>;

        dirNode->setPivot(pivot);
        dirNode->setMu(mu);
        dirNode->setCoverage(Dataset<type>::getMax(*dataset, pivot, df));

        //getNodeHistogram(*dataset, pivot, df, dirNode);

        incrementDistanceVector(pivot, previousPivots.back(), &distanceVector, df);
        previousPivots.push_back(pivot);

        Dataset<type> *innerDataset = new Dataset<type>;
        Dataset<type> *outerDataset = new Dataset<type>;
        splitDataset(pivot, mu, *dataset, innerDataset, outerDataset, df);

        delete(dataset);

        if (innerDataset->getSize() > 0){
            dirNode->setLeftNode(buildTree(true, threshold, maxElements, previousPivots, distanceVector,
                                           innerDataset, df/*, progressBar*/));
        } else {
            delete(innerDataset);
            dirNode->setLeftNode(nullptr);
        }

        if (outerDataset->getSize() > 0){
            dirNode->setRightNode(buildTree(false, threshold, maxElements, previousPivots, distanceVector,
                                            outerDataset, df/*, progressBar*/));
        } else {
            delete(outerDataset);
            dirNode->setRightNode(nullptr);
        }

        return dirNode;

    } else {

        LeafNodeVPTree<BasicArrayObject<type>> *leaf = new LeafNodeVPTree<BasicArrayObject<type>>;

        leaf->setCircunscribed(isCircunscribed);
        leaf->setPivot(pivot);
        leaf->setMu(mu);

        leaf->setCoverage(Dataset<type>::getMax(*dataset, pivot, df));

        incrementDistanceVector(pivot, previousPivots.back(), &distanceVector, df);
        leaf->setDistanceVector(distanceVector);

        previousPivots.push_back(pivot);
        leaf->setPreviousPivots(previousPivots);

//        dataset->erase(pivot);
//        updateProgress(progressBar);

        for (size_t i = 0; i < dataset->getSize(); i++){

            //leaf->push_back(Pair<BasicArrayObject<type>>(dataset->operator [](i), distanceVector.back()));
            leaf->push_back(Pair<BasicArrayObject<type>>(*dataset->instance(i), distanceVector.back()));
            //updateProgress(progressBar);
        }

        delete(dataset);

//        if (dataset->getCardinality() > 0)
//            getBucket(threshold, previousPivots, distanceVector, dataset, df, &leaf->getBucket(),
//                      progressBar);

        return leaf;
    }
}

IMPL_TEMPL void VpTree<type, DistanceFunction>::getBucket(const double threshold,
                                                           std::vector<BasicArrayObject<type>> previousPivots,
                                                           std::vector<std::vector<double>> distanceVector,
                                                           Dataset<type> *dataset,
                                                           DistanceFunction *df,
                                                           Bucket<BasicArrayObject<type>> *bucket/*,
                                                           QProgressBar *progressBar*/){


//    updateProgress(progressBar);

    srand(time(nullptr));
    const BasicArrayObject<type> pivot = dataset->operator [](rand()% dataset->getCardinality());

    const double mu = dataset->medianDistance(pivot, df);

    if ((mu > threshold) && (dataset->getCardinality() > 1)){

        // Split the dataset in innerDataset (< mu) and outerDataset (>= mu).
        Dataset<type> *innerDataset = new Dataset<type>;
        Dataset<type> *outerDataset = new Dataset<type>;
        splitDataset(pivot, mu, *dataset, innerDataset, outerDataset, df);

        // Deletes the dataset to save memory.
        delete(dataset);

        // Increments the DistanceVector with a vector of distances for this node.
        incrementDistanceVector(pivot, previousPivots.back(), &distanceVector, df);

        bucket->push_back(Pair<BasicArrayObject<type>>(pivot, distanceVector.back()));

        previousPivots.push_back(pivot);

        // Recursively calls itself with the inner part of the dataset.
        if (innerDataset->getSize() > 0){
            getBucket(threshold, previousPivots, distanceVector, innerDataset, df, bucket/*, progressBar*/);
        } else {
            delete(innerDataset);
        }

        // Recursively calls itself with the outer part of the dataset.
        if (outerDataset->getSize() > 0){
            getBucket(threshold, previousPivots, distanceVector, outerDataset, df, bucket/*, progressBar*/);
        } else {
            delete(outerDataset);
        }

    } else {
        delete(dataset);

        incrementDistanceVector(pivot, previousPivots.back(), &distanceVector, df);

        bucket->push_back(Pair<BasicArrayObject<type>>(pivot, distanceVector.back()));
    }

}

IMPL_TEMPL void VpTree<type, DistanceFunction>::splitDataset(const BasicArrayObject<type> &pivot,
                                                              const double mu,
                                                              const Dataset<type> &sourceDataset,
                                                              Dataset<type> *innerDataset,
                                                              Dataset<type> *outerDataset,
                                                              DistanceFunction *df){

    const size_t cardinality = sourceDataset.getCardinality();

//    const uint32_t pivotOID = pivot.getOID();

    for (size_t i = 0; i < cardinality; i++){

        const BasicArrayObject<type> &fv = sourceDataset.getFeatureVector(i);

//        if (pivotOID != fv.getOID()){
            if (df->getDistance(pivot, fv) < mu)
                innerDataset->push_back(fv);
            else
                outerDataset->push_back(fv);
//        }

    }

    innerDataset->setCardinality(innerDataset->getSize());
    outerDataset->setCardinality(outerDataset->getSize());

    innerDataset->setDimensionality(sourceDataset.getDimensionality());
    outerDataset->setDimensionality(sourceDataset.getDimensionality());
}

IMPL_TEMPL void VpTree<type, DistanceFunction>::incrementDistanceVector(const BasicArrayObject<type> &pivot, const BasicArrayObject<type> &parentPivot,
                                                                         std::vector<std::vector<double>> *distanceVector,
                                                                         DistanceFunction *df) const{

    std::vector<double> distances;

    if (distanceVector->begin() != distanceVector->end()){

        distances = distanceVector->back();

        const double distanceToParent = df->getDistance(pivot, parentPivot);

        // Increments every element of distances by distanceToParent.
        std::for_each(distances.begin(), distances.end(), [distanceToParent](double &dist){
            dist += distanceToParent;
        });

        distances.push_back(distanceToParent);
    }
    distanceVector->push_back(distances);

}

IMPL_TEMPL size_t VpTree<type, DistanceFunction>::getHeight(Node<BasicArrayObject<type>> *node) const{

    if (node == nullptr) return -1;

    DirectorNode<BasicArrayObject<type>> *dirNode = dynamic_cast<DirectorNode<BasicArrayObject<type>>*>(node);

    if (dirNode != nullptr)
        return 1 + std::max(getHeight(dirNode->getLeftNode()), getHeight(dirNode->getRightNode()));
    else
        return 1;
}

IMPL_TEMPL void VpTree<type, DistanceFunction>::getLeafs(Node<BasicArrayObject<type>> *node, std::vector<LeafNodeVPTree<BasicArrayObject<type>>*>*leafVec) const{

    if (node != nullptr){

        const DirectorNode<BasicArrayObject<type>> *dirNode = dynamic_cast<DirectorNode<BasicArrayObject<type>>*>(node);

        if (dirNode != nullptr){
            getLeafs(dirNode->getLeftNode(), leafVec);
            getLeafs(dirNode->getRightNode(), leafVec);
        } else {
            leafVec->push_back(dynamic_cast<LeafNodeVPTree<BasicArrayObject<type>>*>(node));
        }

    }

}

IMPL_TEMPL void VpTree<type, DistanceFunction>::getDirectorNodes(Node<BasicArrayObject<type>> *node, std::vector<DirectorNode<BasicArrayObject<type>>*> *ids){

    if (node != nullptr){

        DirectorNode<BasicArrayObject<type>> *dirNode = dynamic_cast<DirectorNode<BasicArrayObject<type>>*>(node);

        if (dirNode != nullptr){

            ids->push_back(dirNode);

            getDirectorNodes(dirNode->getLeftNode(), ids);
            getDirectorNodes(dirNode->getRightNode(), ids);
        }

    }

}

IMPL_TEMPL void VpTree<type, DistanceFunction>::getDirectorNodesByHeight(Node<BasicArrayObject<type>> *node,
        std::vector<std::vector<DirectorNode<BasicArrayObject<type>>*>> *dirVec,
        size_t height) const{

    if (node != nullptr){

        DirectorNode<BasicArrayObject<type>> *dirNode = dynamic_cast<DirectorNode<BasicArrayObject<type>>*>(node);

        if (dirNode != nullptr){

            dirVec->at(height).push_back(dirNode);

            getDirectorNodesByHeight(dirNode->getLeftNode(), dirVec, height+1);
            getDirectorNodesByHeight(dirNode->getRightNode(), dirVec, height+1);
        }

    }
}

//IMPL_TEMPL void VpTree<type, DistanceFunction>::getNodeHistogram(const Dataset<type> &dataset,
//                                                                  const BasicArrayObject<DType> &pivot,
//                                                                  DistanceFunction *df,
//                                                                  DirectorNode<BasicArrayObject<DType>>* dirNode) const{

//    // 20 é o número de bins do histograma.
//    dirNode->setHistogram(PivotHistogram(dataset, pivot, 20, df));
//}

IMPL_TEMPL LeafNodeVPTree<BasicArrayObject<type>> *VpTree<type, DistanceFunction>::generateBalancedLeaf(size_t i, size_t j, std::vector<std::vector<DirectorNode<BasicArrayObject<type>>*>> heightList){

    LeafNodeVPTree<BasicArrayObject<type>> *newLeaf = new LeafNodeVPTree<BasicArrayObject<type>>;

    newLeaf->setPivot(heightList[i+1][j]->getPivot());
    newLeaf->setMu(heightList[i+1][j]->getMu());
    newLeaf->setCircunscribed(true);

    pushExceededNodes(heightList[i+1][j]->getLeftNode(), newLeaf);
    pushExceededNodes(heightList[i+1][j]->getRightNode(), newLeaf);

    return newLeaf;
}

IMPL_TEMPL void VpTree<type, DistanceFunction>::balanceTree(){

    int height = getHeight(root);
    std::vector<std::vector<DirectorNode<BasicArrayObject<type>>*>> heightList(height);
    getDirectorNodesByHeight(root, &heightList);

    for (size_t i = 0; i < heightList.size()-1; i++){

        // If next height is less than double the size of actual height, the next height is unbalanced.
        if (heightList[i].size()*2 > heightList[i+1].size()){

            // For every node in the unbalanced height.
            for (size_t j = 0; j < heightList[i+1].size(); j++){

                // Search for nodes pointing to the unbalanced node and points it to the new node.
                for (size_t k = 0; k < heightList[i].size(); k++){

                    if (heightList[i][k]->getLeftNode() == heightList[i+1][j])
                        heightList[i][k]->setLeftNode(generateBalancedLeaf(i, j, heightList));

                    if (heightList[i][k]->getRightNode() == heightList[i+1][j])
                        heightList[i][k]->setRightNode(generateBalancedLeaf(i, j, heightList));
                }
            }

            break;
        }
    }
}

IMPL_TEMPL void VpTree<type, DistanceFunction>::pushExceededNodes(Node<BasicArrayObject<type>> *node,
                                                                   LeafNodeVPTree<BasicArrayObject<type>> *newLeaf){

    if (node == nullptr) return;

    newLeaf->push_exceeded_node(node);

    DirectorNode<BasicArrayObject<type>> *dirNode = dynamic_cast<DirectorNode<BasicArrayObject<type>>*>(node);

    if (dirNode != nullptr){

        pushExceededNodes(dirNode->getLeftNode(), newLeaf);
        pushExceededNodes(dirNode->getRightNode(), newLeaf);

    } else {

        LeafNodeVPTree<BasicArrayObject<type>> *leafNode = static_cast<LeafNodeVPTree<BasicArrayObject<type>>*>(node);

        for (size_t i = 0; i < leafNode->getBucket().numberOfElements(); i++){
            newLeaf->push_back(leafNode->getPair(i));
        }
    }
}

/**
 *@brief A getter to the root.
 *@returns a pointer to the root of the tree
 */
IMPL_TEMPL Node<BasicArrayObject<type>> *VpTree<type, DistanceFunction>::getRoot() const{ return root; }

IMPL_TEMPL void VpTree<type, DistanceFunction>::rangeQuery(const BasicArrayObject<type> &qElement,
                                                            double radius,
                                                            Node<BasicArrayObject<type>> *node,
                                                            Dataset<type> *answer,
                                                            DistanceFunction *df){


     if (node == nullptr) return;
     node->wasVisited = true;

//JOAO - CODE
//     if (df->getDistance(qElement, node->getPivot()) <= radius)
//         answer->push_back(node->getPivot());

     LeafNodeVPTree<BasicArrayObject<type>> *leaf = dynamic_cast<LeafNodeVPTree<BasicArrayObject<type>>*>(node);

     if (leaf != nullptr) {

         for (size_t i = 0; i < leaf->numberOfElements(); i++){

             if (df->getDistance(qElement, leaf->getPair(i).first()) <= radius)
                 answer->push_back(leaf->getPair(i).first());
         }

     } else {

        double dist = df->getDistance(node->getPivot(), qElement);

        DirectorNode<BasicArrayObject<type>> *dirNode = dynamic_cast<DirectorNode<BasicArrayObject<type>>*>(node);

        if (dist < node->getMu()) {

            if ((dist - radius) <= node->getMu())
                rangeQuery(qElement, radius, dirNode->getLeftNode(), answer, df);

            if ((dist + radius) >= node->getMu())
                rangeQuery(qElement, radius, dirNode->getRightNode(), answer, df);

        } else {

            if ((dist + radius) >= node->getMu())
                rangeQuery(qElement, radius, dirNode->getRightNode(), answer, df);

            if ((dist - radius) <= node->getMu())
                rangeQuery(qElement, radius, dirNode->getLeftNode(), answer, df);
        }

     }

}

enum PositionRelativeToPartition{INSIDE, RING, OUTSIDE};

IMPL_TEMPL auto VpTree<type, DistanceFunction>::checkSqPosition(double dist, double mu, double max) {
    if (dist < mu)
        return PositionRelativeToPartition::INSIDE;
    else if (dist <= max)
        return PositionRelativeToPartition::RING;
    else
        return PositionRelativeToPartition::OUTSIDE;
}

IMPL_TEMPL void VpTree<type, DistanceFunction>::setLeftRightMinMax(const BasicArrayObject<type> &qElement,
                                                                    DirectorNode<BasicArrayObject<type>> *dirNode,
                                                                    DistanceFunction *df){


    double dist =  df->getDistance(qElement, dirNode->getPivot());

    PositionRelativeToPartition sqPosition = checkSqPosition(dist, dirNode->getMu(), dirNode->getCoverage());

    // sq dentro bola
    if (sqPosition == PositionRelativeToPartition::INSIDE) {

        // Min do Left
        dirNode->getLeftNode()->setDMin(0.0);

        // Max do Left
        if (dirNode->getDMax() > dirNode->getMu() + dist)
            dirNode->getLeftNode()->setDMax(dirNode->getMu() + dist);
        else
            dirNode->getLeftNode()->setDMax(dirNode->getDMax());

        // Min do Right
        dirNode->getRightNode()->setDMin(dirNode->getMu() - dist);

        // Max do Right
        if (dirNode->getDMax() > dirNode->getCoverage() + dist)
            dirNode->getRightNode()->setDMax(dirNode->getCoverage() + dist);
        else
            dirNode->getRightNode()->setDMax(dirNode->getDMax());
    }
    // sq no anel
    else if (sqPosition == PositionRelativeToPartition::RING) {

        // Min do Left
        dirNode->getLeftNode()->setDMin(dist - dirNode->getMu());

        // Max do Left
        if (dirNode->getDMax() > dist + dirNode->getMu())
            dirNode->getLeftNode()->setDMax(dist + dirNode->getMu());
        else
            dirNode->getLeftNode()->setDMax(dirNode->getDMax());

        // Min do Right
        dirNode->getRightNode()->setDMin(0.0);

        // Max do Right
        if (dirNode->getDMax() > dist + dirNode->getCoverage())
            dirNode->getRightNode()->setDMax(dist + dirNode->getCoverage());
        else
            dirNode->getRightNode()->setDMax(dirNode->getDMax());
    }
    // sq fora da partição
    else {

        // Min do Left
        dirNode->getLeftNode()->setDMin(dist - dirNode->getMu());

        // Max do Left
        if (dirNode->getDMax() > dist + dirNode->getMu())
            dirNode->getLeftNode()->setDMax(dist + dirNode->getMu());
        else
            dirNode->getLeftNode()->setDMax(dirNode->getDMax());

        // Min do Right
        dirNode->getRightNode()->setDMin(dist - dirNode->getCoverage());

        // Max do Right
        if (dirNode->getDMax() > dist + dirNode->getCoverage())
            dirNode->getRightNode()->setDMax(dist + dirNode->getCoverage());
        else
            dirNode->getRightNode()->setDMax(dirNode->getDMax());
    }
}

template <class T>
class NodeCmpFromPtr{

public:
    bool operator() (Node<T> *node1, Node<T> *node2){
        return node1->operator>(*node2);
    }
};



IMPL_TEMPL void VpTree<type, DistanceFunction>::kNNIncWalk(const BasicArrayObject<type> &qElement,
                                                            size_t k,
                                                            Node<BasicArrayObject<type>> *tree,
                                                            std::priority_queue<QueueItem<BasicArrayObject<type>>,
                                                                                std::vector<QueueItem<BasicArrayObject<type>>>,
                                                                                std::greater<QueueItem<BasicArrayObject<type>>>> *resultQueue,
                                                            DistanceFunction *df){

    std::priority_queue<Node<BasicArrayObject<type>> *,
                        std::vector<Node<BasicArrayObject<type>> *>,
                        NodeCmpFromPtr<BasicArrayObject<type>>> nodeQueue;

    std::priority_queue<QueueItem<BasicArrayObject<type>>,
                        std::vector<QueueItem<BasicArrayObject<type>>>,
                        std::greater<QueueItem<BasicArrayObject<type>>> > elementQueue;

    tree->setDMin(0.0);
    tree->setDMax(std::numeric_limits<double>::infinity());

    nodeQueue.push(tree);

    Dataset<type>* datasetLeaf;

    // Não completou resultado && há nós ou elementos a serem obtidos.
    while (resultQueue->size() < k && !(nodeQueue.empty() && elementQueue.empty())){

        // Se entrar aqui, é garantido que o nodeQueue tem nós e que o elementQueue está vazio.
        if (elementQueue.empty()) {

            // Popa a partição
            Node<BasicArrayObject<type>> *node = nodeQueue.top();
            nodeQueue.pop();

            // Se o cast falhar, nó não é do tipo Leaf e retorna nullptr.
            LeafNodeVPTree<BasicArrayObject<type>> *leaf = dynamic_cast<LeafNodeVPTree<BasicArrayObject<type>> *>(node);
            if (leaf != nullptr) {

                leafNodeAccess++;

                datasetLeaf = readLeafNode(node);

                for(size_t i = 0; i < datasetLeaf->getCardinality(); i++)
                {

                    QueueItem<BasicArrayObject<type>> qi = QueueItem<BasicArrayObject<type>>(df->getDistance(qElement, datasetLeaf->getFeatureVector(i)), datasetLeaf->getFeatureVector(i));
                    elementQueue.push(qi);

                }

                delete datasetLeaf;

                // Insere os elementos do folha no elementQueue
//                for (int i = 0; i < leaf->numberOfElements(); i++){
//                    BasicArrayObject<type> leafElement = leaf->getPair(i).first();
//                    double dist = df->getDistance(leafElement, qElement);
//                    QueueItem<BasicArrayObject<type>> qi = QueueItem<BasicArrayObject<type>>(dist, leafElement);
//                    elementQueue.push(qi);
//                }

//                // Seta variável para verificarmos quais partições foram visitadas.
                node->wasVisited = true;

            } else {

                DirectorNode<BasicArrayObject<type>> *dirNode = dynamic_cast<DirectorNode<BasicArrayObject<type>> *>(node);

                // Seta min/max dos filhos
                setLeftRightMinMax(qElement, dirNode, df);

                // Insere os filhos no nodeQueue
                nodeQueue.push(dirNode->getLeftNode());
                nodeQueue.push(dirNode->getRightNode());

            }

            // Se entrar aqui, é garantido que o elementQueue tem elementos e que o nodeQueue tem nós,
            // assim como que o dMin do topo do nodeQueue é menor do que a distância no topo do elementQueue.
        } else if (!nodeQueue.empty() &&
                   nodeQueue.top()->getDMin() < elementQueue.top().dist) {

            // Popa a partição
            Node<BasicArrayObject<type>> *node = nodeQueue.top();
            nodeQueue.pop();

            // Checa se é folha
            // Se o cast falhar, nó não é do tipo Leaf e retorna nullptr.
            LeafNodeVPTree<BasicArrayObject<type>> *leaf = dynamic_cast<LeafNodeVPTree<BasicArrayObject<type>>*>(node);
            if (leaf != nullptr) {

                leafNodeAccess++;

                datasetLeaf = readLeafNode(node);

                for(size_t i = 0; i < datasetLeaf->getCardinality(); i++)
                {

                    QueueItem<BasicArrayObject<type>> qi = QueueItem<BasicArrayObject<type>>(df->getDistance(qElement, datasetLeaf->getFeatureVector(i)), datasetLeaf->getFeatureVector(i));
                    elementQueue.push(qi);

                }

                delete datasetLeaf;

                // Insere os elementos do folha no elementQueue
//                for (size_t i = 0; i < leaf->numberOfElements(); i++){
//                    BasicArrayObject<type> leafElement = leaf->getPair(i).first();
//                    double dist = df->getDistance(leafElement, qElement);
//                    QueueItem<BasicArrayObject<type>> qi = QueueItem<BasicArrayObject<type>>(dist, leafElement);
//                    elementQueue.push(qi);
//                }

                // Seta variável para verificarmos quais partições foram visitadas.
                node->wasVisited = true;

            } else {

                DirectorNode<BasicArrayObject<type>> *dirNode = dynamic_cast<DirectorNode<BasicArrayObject<type>>*>(node);

                // Seta min/max dos filhos
                setLeftRightMinMax(qElement, dirNode, df);

                // Insere os filhos no nodeQueue
                nodeQueue.push(dirNode->getLeftNode());
                nodeQueue.push(dirNode->getRightNode());

            }

            // Se entrar aqui, é garantido ou que __nodeQueue está vazio, mas o elementQueue não__,
            // ou que __a distância do elemento no topo do elementQueue é menor do que o dMin do topo do nodeQueue__
        } else {
            resultQueue->push(elementQueue.top());
            elementQueue.pop();
        }
    }
}

IMPL_TEMPL void VpTree<type, DistanceFunction>::kNNInc(const BasicArrayObject<type> &qElement,
                                                     size_t k,
                                                     Node<BasicArrayObject<type>> *node,
                                                     Dataset<type> *answer,
                                                     DistanceFunction *df){

    std::priority_queue<QueueItem<BasicArrayObject<type>>,
                        std::vector<QueueItem<BasicArrayObject<type>>>,
                        std::greater<QueueItem<BasicArrayObject<type>>>> queue;

    leafNodeAccess = 0;
    df->resetStatistics();

    kNNIncWalk(qElement, k, node, &queue, df);

    //std::cout << "(kNNInc) resultQueue.size(): " << queue.size() << std::endl;
    answer->reserve(k);
    for (size_t i = 0; i < k; i++){
        answer->push_back(queue.top().featureVector);
        //std::cout << "(kNNInc) Distancia do " << i << " vizinho mais proximo:" << queue.top().dist << std::endl;
        queue.pop();
    }
}


IMPL_TEMPL void VpTree<type, DistanceFunction>::kNNWalk(const BasicArrayObject<type> &qElement,
                                                         double &tau,
                                                         size_t k,
                                                         Node<BasicArrayObject<type>> *node,
                                                         std::priority_queue<QueueItem<BasicArrayObject<type>>,
                                                                             std::vector<QueueItem<BasicArrayObject<type>>>,
                                                                             std::less<QueueItem<BasicArrayObject<type>>>> *queue,
                                                         std::vector<std::pair<int, double>> *pivotVec,
                                                         DistanceFunction *df){

    node->wasVisited = true;

    LeafNodeVPTree<BasicArrayObject<type>> *leaf = dynamic_cast<LeafNodeVPTree<BasicArrayObject<type>>*>(node);

    if (leaf != nullptr) {

        for (size_t i = 0; i < leaf->numberOfElements(); i++){

            BasicArrayObject<type> leafElement = leaf->getPair(i).first();
            double dist = 0.0;

            bool found = 0;
            for (size_t j = 0; j < pivotVec->size(); j++){
                if (pivotVec->operator[](j).first == leafElement.getOID()){
                    dist = pivotVec->operator[](j).second;
                    found = 1;
                    break;
                }
            }

            if (!found)
                dist = df->getDistance(qElement, leafElement);

            queue->push(QueueItem<BasicArrayObject<type>>(dist, leafElement));

            if (queue->size() > k)
                queue->pop();
        }

    } else {

        DirectorNode<BasicArrayObject<type>> *dirNode = dynamic_cast<DirectorNode<BasicArrayObject<type>>*>(node);

        double dist = df->getDistance(dirNode->getPivot(), qElement);

        pivotVec->push_back(std::pair<int, double>(dirNode->getPivot().getOID(), dist));

//        if (dist <= tau){

//            queue->push(QueueItem<BasicArrayObject<DType>>(dist, dirNode->getPivot()));

//            while (queue->size() > k)
//                queue->pop();

//            if (queue->size() == k)
//                tau = queue->top().dist;
//        }

        if (dist < node->getMu()){

            if ((dist - tau <= dirNode->getMu()) && dirNode->getLeftNode() != nullptr)
                kNNWalk(qElement, tau, k, dirNode->getLeftNode(), queue, pivotVec, df);

            if ((dist + tau >= dirNode->getMu()) && dirNode->getRightNode() != nullptr)
                kNNWalk(qElement, tau, k, dirNode->getRightNode(), queue, pivotVec, df);

        } else {

            if ((dist + tau >= dirNode->getMu()) && dirNode->getRightNode() != nullptr)
                kNNWalk(qElement, tau, k, dirNode->getRightNode(), queue, pivotVec, df);

            if ((dist - tau <= node->getMu()) && dirNode->getLeftNode() != nullptr)
                kNNWalk(qElement, tau, k, dirNode->getLeftNode(), queue, pivotVec, df);

        }

    }
}

IMPL_TEMPL void VpTree<type, DistanceFunction>::kNN(const BasicArrayObject<type> &qElement,
                                                     size_t k,
                                                     Node<BasicArrayObject<type>> *node,
                                                     Dataset<type> *answer,
                                                     DistanceFunction *df){

    std::priority_queue<QueueItem<BasicArrayObject<type>>,
                        std::vector<QueueItem<BasicArrayObject<type>>>,
                        std::less<QueueItem<BasicArrayObject<type>>>> queue;


    double tau = std::numeric_limits<double>::infinity();
    kNNWalk(qElement, tau, k, node, &queue, new std::vector<std::pair<int, double>>(), df);
    //JOAO - CODE
    //kNNWalk(qElement, tau, k, node, &queue, df);


    std::priority_queue<QueueItem<BasicArrayObject<type>>,
                        std::vector<QueueItem<BasicArrayObject<type>>>,
                        std::greater<QueueItem<BasicArrayObject<type>>>> reverse;

    for (size_t i = 0; i < k; i++){
        reverse.push(queue.top());
        queue.pop();
    }

    //std::cout << "(kNN Padrão) resultQueue.size(): " << reverse.size() << std::endl;

    answer->reserve(k);
    for (size_t i = 0; i < k; i++){
//        answer->push_back(reverse.top().featureVector);
        answer[k-i] = reverse.top().featureVector;
        //std::cout << "(kNN Padrao) Distancia do " << i << " vizinho mais proximo:" << reverse.top().dist << std::endl;
        reverse.pop();
    }


}

//IMPL_TEMPL void VpTree<BasicArrayObject<DType>, DistanceFunction>::updateProgress(QProgressBar *progressBar){
////    ++progress;
////    if (progressBar->maximum() < SHRT_MAX)
////        emit progressBarValueChanged(progress);
////    else if (progress % SHRT_MAX == 0)
////        emit progressBarValueChanged(progress/SHRT_MAX);
//}
































#endif // VPTREE_H
