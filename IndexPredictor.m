function out = IndexPredictor(features,models)
Index = ["Wemmert Gancarski","Point Biserial","AUCC","Silhouette","CDbw","VRC","DB","WB","DBCV"];
nalg = size(models.knn,2);
prediction = zeros(1,nalg);

for i = 1:nalg
    prediction(i) = predict(models.knn{i},features');
end

if sum(prediction) > 1
    [~,sel] = max(bsxfun(@times,prediction,models.precision'),[],2);
end

out.goodIndex = Index(prediction==1);
out.selectedIndex = Index(sel);

end