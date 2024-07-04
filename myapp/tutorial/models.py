from django.db import models

# Create your models here.
class Cell(models.Model):
    name = models.CharField(max_length=100)

    def __str__(self):
        return f"This is cell instance of {self.name}"

class Gene(models.Model):
    name = models.CharField(max_length=100)

    def __str__(self):
        return f"This is gene instance of {self.name}"
